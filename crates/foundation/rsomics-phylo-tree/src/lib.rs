#![allow(clippy::missing_errors_doc, clippy::missing_panics_doc)]

use std::fmt::Write as _;

pub type NodeId = usize;

#[derive(Debug, Clone)]
pub struct Node {
    pub id: NodeId,
    pub name: Option<String>,
    pub branch_length: Option<f64>,
    pub parent: Option<NodeId>,
    pub children: Vec<NodeId>,
}

#[derive(Debug, Clone, Default)]
pub struct Tree {
    pub nodes: Vec<Node>,
    pub root: NodeId,
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum NewickError {
    #[error("unexpected character {ch:?} at offset {pos}")]
    Unexpected { ch: char, pos: usize },
    #[error("unbalanced parentheses (depth {depth} at end)")]
    Unbalanced { depth: i64 },
    #[error("malformed branch length {token:?}: {reason}")]
    BadBranchLength { token: String, reason: String },
    #[error("empty input")]
    Empty,
}

pub type Result<T> = std::result::Result<T, NewickError>;

impl Tree {
    pub fn from_newick(src: &str) -> Result<Self> {
        let trimmed = src.trim();
        if trimmed.is_empty() {
            return Err(NewickError::Empty);
        }
        let mut parser = Parser {
            src: trimmed.as_bytes(),
            pos: 0,
        };
        let mut tree = Tree::default();
        let root = parser.parse_subtree(&mut tree, None)?;
        tree.root = root;
        // Newick trees terminate with `;`. Tolerate missing-semicolon to be
        // generous; surface anything past the terminator as an error.
        parser.skip_ws();
        if parser.pos < parser.src.len() && parser.src[parser.pos] == b';' {
            parser.pos += 1;
        }
        parser.skip_ws();
        if parser.pos != parser.src.len() {
            return Err(NewickError::Unexpected {
                ch: parser.src[parser.pos] as char,
                pos: parser.pos,
            });
        }
        Ok(tree)
    }

    #[must_use]
    pub fn to_newick(&self) -> String {
        let mut out = String::new();
        self.write_node(self.root, &mut out);
        out.push(';');
        out
    }

    fn write_node(&self, id: NodeId, out: &mut String) {
        let n = &self.nodes[id];
        if !n.children.is_empty() {
            out.push('(');
            for (i, &c) in n.children.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                self.write_node(c, out);
            }
            out.push(')');
        }
        if let Some(name) = &n.name {
            out.push_str(name);
        }
        if let Some(bl) = n.branch_length {
            let _ = write!(out, ":{bl}");
        }
    }

    pub fn leaves(&self) -> impl Iterator<Item = &Node> {
        self.nodes.iter().filter(|n| n.children.is_empty())
    }

    #[must_use]
    pub fn n_leaves(&self) -> usize {
        self.leaves().count()
    }
}

struct Parser<'a> {
    src: &'a [u8],
    pos: usize,
}

impl Parser<'_> {
    fn skip_ws(&mut self) {
        while self.pos < self.src.len() && self.src[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }

    fn parse_subtree(&mut self, tree: &mut Tree, parent: Option<NodeId>) -> Result<NodeId> {
        self.skip_ws();
        let id = tree.nodes.len();
        tree.nodes.push(Node {
            id,
            name: None,
            branch_length: None,
            parent,
            children: Vec::new(),
        });

        if self.pos < self.src.len() && self.src[self.pos] == b'(' {
            self.pos += 1;
            loop {
                let child = self.parse_subtree(tree, Some(id))?;
                tree.nodes[id].children.push(child);
                self.skip_ws();
                match self.src.get(self.pos) {
                    Some(&b',') => {
                        self.pos += 1;
                    }
                    Some(&b')') => {
                        self.pos += 1;
                        break;
                    }
                    Some(&ch) => {
                        return Err(NewickError::Unexpected {
                            ch: ch as char,
                            pos: self.pos,
                        });
                    }
                    None => return Err(NewickError::Unbalanced { depth: 1 }),
                }
            }
        }

        // Optional name + branch length.
        let name = self.read_token();
        if !name.is_empty() {
            tree.nodes[id].name = Some(name);
        }
        self.skip_ws();
        if self.pos < self.src.len() && self.src[self.pos] == b':' {
            self.pos += 1;
            let bl = self.read_token();
            let v: f64 = bl.parse().map_err(|e: std::num::ParseFloatError| {
                NewickError::BadBranchLength {
                    token: bl.clone(),
                    reason: e.to_string(),
                }
            })?;
            tree.nodes[id].branch_length = Some(v);
        }
        Ok(id)
    }

    fn read_token(&mut self) -> String {
        let start = self.pos;
        while let Some(&ch) = self.src.get(self.pos) {
            if ch == b'('
                || ch == b')'
                || ch == b','
                || ch == b':'
                || ch == b';'
                || ch.is_ascii_whitespace()
            {
                break;
            }
            self.pos += 1;
        }
        std::str::from_utf8(&self.src[start..self.pos])
            .unwrap_or("")
            .to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_balanced_tree() {
        // 4-leaf bifurcating: ((A,B),(C,D));
        let tree = Tree::from_newick("((A,B),(C,D));").unwrap();
        assert_eq!(tree.n_leaves(), 4);
        let names: Vec<&str> = tree
            .leaves()
            .map(|n| n.name.as_deref().unwrap_or(""))
            .collect();
        for want in ["A", "B", "C", "D"] {
            assert!(names.contains(&want), "missing leaf {want}: {names:?}");
        }
    }

    #[test]
    fn parse_branch_lengths() {
        let tree = Tree::from_newick("((A:0.1,B:0.2):0.3,C:0.4);").unwrap();
        let a = tree
            .leaves()
            .find(|n| n.name.as_deref() == Some("A"))
            .unwrap();
        let c = tree
            .leaves()
            .find(|n| n.name.as_deref() == Some("C"))
            .unwrap();
        assert!((a.branch_length.unwrap() - 0.1).abs() < 1e-9);
        assert!((c.branch_length.unwrap() - 0.4).abs() < 1e-9);
    }

    #[test]
    fn round_trip_simple() {
        let src = "((A:0.1,B:0.2),(C,D:0.4));";
        let tree = Tree::from_newick(src).unwrap();
        let out = tree.to_newick();
        let again = Tree::from_newick(&out).unwrap();
        assert_eq!(again.n_leaves(), tree.n_leaves());
    }

    #[test]
    fn parse_single_leaf() {
        let tree = Tree::from_newick("A;").unwrap();
        assert_eq!(tree.n_leaves(), 1);
        assert_eq!(tree.nodes[tree.root].name.as_deref(), Some("A"));
    }

    #[test]
    fn parse_missing_terminator_tolerated() {
        let tree = Tree::from_newick("(A,B)").unwrap();
        assert_eq!(tree.n_leaves(), 2);
    }

    #[test]
    fn empty_input_rejected() {
        assert!(matches!(Tree::from_newick("   "), Err(NewickError::Empty)));
    }

    #[test]
    fn unbalanced_paren_rejected() {
        assert!(Tree::from_newick("((A,B);").is_err());
    }

    #[test]
    fn bad_branch_length_rejected() {
        assert!(matches!(
            Tree::from_newick("(A:notanumber,B);"),
            Err(NewickError::BadBranchLength { .. })
        ));
    }
}
