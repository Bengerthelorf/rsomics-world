#![allow(clippy::missing_errors_doc, clippy::missing_panics_doc)]

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub type TaxId = u32;

#[derive(Debug, Clone)]
pub struct TaxNode {
    pub tax_id: TaxId,
    pub parent: TaxId,
    pub rank: String,
    pub name: String,
}

#[derive(Debug, Default, Clone)]
pub struct Taxonomy {
    pub nodes: HashMap<TaxId, TaxNode>,
    pub merged: HashMap<TaxId, TaxId>,
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum TaxonomyError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("malformed taxdump line {line} in {file}: {reason}")]
    Malformed {
        file: String,
        line: usize,
        reason: String,
    },
    #[error("unknown tax_id: {0}")]
    Unknown(TaxId),
}

pub type Result<T> = std::result::Result<T, TaxonomyError>;

impl Taxonomy {
    pub fn from_dump(nodes_dmp: &Path, names_dmp: &Path) -> Result<Self> {
        let mut nodes: HashMap<TaxId, TaxNode> = HashMap::new();
        for (lineno, line) in BufReader::new(File::open(nodes_dmp)?).lines().enumerate() {
            let line = line?;
            let fields = parse_dump_line(&line);
            if fields.len() < 3 {
                return Err(TaxonomyError::Malformed {
                    file: nodes_dmp.display().to_string(),
                    line: lineno + 1,
                    reason: "fewer than 3 fields".into(),
                });
            }
            let tax_id: TaxId = fields[0].parse().map_err(|_| TaxonomyError::Malformed {
                file: nodes_dmp.display().to_string(),
                line: lineno + 1,
                reason: format!("bad tax_id {:?}", fields[0]),
            })?;
            let parent: TaxId = fields[1].parse().map_err(|_| TaxonomyError::Malformed {
                file: nodes_dmp.display().to_string(),
                line: lineno + 1,
                reason: format!("bad parent_tax_id {:?}", fields[1]),
            })?;
            nodes.insert(
                tax_id,
                TaxNode {
                    tax_id,
                    parent,
                    rank: fields[2].to_string(),
                    name: String::new(),
                },
            );
        }
        for (lineno, line) in BufReader::new(File::open(names_dmp)?).lines().enumerate() {
            let line = line?;
            let fields = parse_dump_line(&line);
            if fields.len() < 4 {
                return Err(TaxonomyError::Malformed {
                    file: names_dmp.display().to_string(),
                    line: lineno + 1,
                    reason: "fewer than 4 fields".into(),
                });
            }
            if fields[3] != "scientific name" {
                // names.dmp also carries synonym/common; skip those
                continue;
            }
            let tax_id: TaxId = fields[0].parse().map_err(|_| TaxonomyError::Malformed {
                file: names_dmp.display().to_string(),
                line: lineno + 1,
                reason: format!("bad tax_id {:?}", fields[0]),
            })?;
            if let Some(node) = nodes.get_mut(&tax_id) {
                node.name = fields[1].to_string();
            }
        }
        Ok(Self {
            nodes,
            merged: HashMap::new(),
        })
    }

    pub fn with_merged(mut self, merged_dmp: &Path) -> Result<Self> {
        for (lineno, line) in BufReader::new(File::open(merged_dmp)?).lines().enumerate() {
            let line = line?;
            let fields = parse_dump_line(&line);
            if fields.len() < 2 {
                return Err(TaxonomyError::Malformed {
                    file: merged_dmp.display().to_string(),
                    line: lineno + 1,
                    reason: "fewer than 2 fields".into(),
                });
            }
            let old: TaxId = fields[0].parse().map_err(|_| TaxonomyError::Malformed {
                file: merged_dmp.display().to_string(),
                line: lineno + 1,
                reason: format!("bad old_tax_id {:?}", fields[0]),
            })?;
            let new: TaxId = fields[1].parse().map_err(|_| TaxonomyError::Malformed {
                file: merged_dmp.display().to_string(),
                line: lineno + 1,
                reason: format!("bad new_tax_id {:?}", fields[1]),
            })?;
            self.merged.insert(old, new);
        }
        Ok(self)
    }

    /// Walk merged-redirects until a live node is reached; returns `id` when
    /// no merge entry exists.
    #[must_use]
    pub fn resolve(&self, id: TaxId) -> TaxId {
        let mut cur = id;
        let mut hops = 0;
        while hops < 16 {
            // cap against circular merged.dmp
            match self.merged.get(&cur) {
                Some(&next) if next != cur => {
                    cur = next;
                    hops += 1;
                }
                _ => break,
            }
        }
        cur
    }

    /// Path from `id` to root (id first, root last).
    pub fn lineage(&self, id: TaxId) -> Result<Vec<TaxId>> {
        let resolved = self.resolve(id);
        if !self.nodes.contains_key(&resolved) {
            return Err(TaxonomyError::Unknown(id));
        }
        let mut path = Vec::new();
        let mut cur = resolved;
        let mut hops = 0;
        loop {
            path.push(cur);
            let node = self.nodes.get(&cur).ok_or(TaxonomyError::Unknown(cur))?;
            if node.parent == cur || hops > 100_000 {
                break;
            }
            cur = node.parent;
            hops += 1;
        }
        Ok(path)
    }

    /// Lowest common ancestor of two or more `tax_ids`.
    pub fn lca(&self, ids: &[TaxId]) -> Result<TaxId> {
        if ids.is_empty() {
            return Err(TaxonomyError::Unknown(0));
        }
        let mut paths: Vec<Vec<TaxId>> = Vec::with_capacity(ids.len());
        for &id in ids {
            paths.push(self.lineage(id)?);
        }
        for p in &mut paths {
            // reverse so index 0 = root, then walk forward
            p.reverse();
        }
        let mut idx = 0;
        loop {
            let candidate = paths[0].get(idx).copied();
            let all_match =
                candidate.is_some_and(|c| paths.iter().all(|p| p.get(idx).copied() == Some(c)));
            if all_match {
                idx += 1;
            } else {
                break;
            }
        }
        let lca_idx = idx.checked_sub(1).ok_or(TaxonomyError::Unknown(0))?;
        Ok(paths[0][lca_idx])
    }
}

fn parse_dump_line(line: &str) -> Vec<&str> {
    // NCBI taxdump: fields separated by `\t|\t`, trailing `\t|`.
    let trimmed = line.strip_suffix("\t|").unwrap_or(line);
    trimmed.split("\t|\t").collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn tax_fixture() -> (tempfile::TempDir, std::path::PathBuf, std::path::PathBuf) {
        let dir = tempfile::tempdir().unwrap();
        let nodes_path = dir.path().join("nodes.dmp");
        let names_path = dir.path().join("names.dmp");
        let mut nodes = std::fs::File::create(&nodes_path).unwrap();
        writeln!(nodes, "1\t|\t1\t|\tno rank\t|").unwrap();
        writeln!(nodes, "2\t|\t1\t|\tsuperkingdom\t|").unwrap();
        writeln!(nodes, "2157\t|\t1\t|\tsuperkingdom\t|").unwrap();
        writeln!(nodes, "1224\t|\t2\t|\tphylum\t|").unwrap();
        writeln!(nodes, "28211\t|\t1224\t|\tclass\t|").unwrap();
        writeln!(nodes, "1236\t|\t1224\t|\tclass\t|").unwrap();
        writeln!(nodes, "28890\t|\t2157\t|\tphylum\t|").unwrap();

        let mut names = std::fs::File::create(&names_path).unwrap();
        writeln!(names, "1\t|\troot\t|\t\t|\tscientific name\t|").unwrap();
        writeln!(names, "2\t|\tBacteria\t|\t\t|\tscientific name\t|").unwrap();
        writeln!(names, "2\t|\tEubacteria\t|\t\t|\tcommon name\t|").unwrap();
        writeln!(names, "2157\t|\tArchaea\t|\t\t|\tscientific name\t|").unwrap();
        writeln!(names, "1224\t|\tProteobacteria\t|\t\t|\tscientific name\t|").unwrap();
        writeln!(
            names,
            "28211\t|\tAlphaproteobacteria\t|\t\t|\tscientific name\t|"
        )
        .unwrap();
        writeln!(
            names,
            "1236\t|\tGammaproteobacteria\t|\t\t|\tscientific name\t|"
        )
        .unwrap();
        writeln!(names, "28890\t|\tEuryarchaeota\t|\t\t|\tscientific name\t|").unwrap();
        (dir, nodes_path, names_path)
    }

    #[test]
    fn from_dump_loads_nodes_and_names() {
        let (_d, nodes, names) = tax_fixture();
        let tax = Taxonomy::from_dump(&nodes, &names).unwrap();
        assert_eq!(tax.nodes.len(), 7);
        assert_eq!(tax.nodes.get(&1224).unwrap().name, "Proteobacteria");
        assert_eq!(tax.nodes.get(&2).unwrap().name, "Bacteria");
        assert_eq!(tax.nodes.get(&1224).unwrap().parent, 2);
    }

    #[test]
    fn lineage_root_to_leaf() {
        let (_d, nodes, names) = tax_fixture();
        let tax = Taxonomy::from_dump(&nodes, &names).unwrap();
        let path = tax.lineage(28211).unwrap();
        assert_eq!(path, vec![28211, 1224, 2, 1]);
    }

    #[test]
    fn lca_within_phylum() {
        let (_d, nodes, names) = tax_fixture();
        let tax = Taxonomy::from_dump(&nodes, &names).unwrap();
        assert_eq!(tax.lca(&[28211, 1236]).unwrap(), 1224);
    }

    #[test]
    fn lca_across_kingdoms_is_root() {
        let (_d, nodes, names) = tax_fixture();
        let tax = Taxonomy::from_dump(&nodes, &names).unwrap();
        assert_eq!(tax.lca(&[28211, 28890]).unwrap(), 1);
    }

    #[test]
    fn merged_redirects_on_resolve() {
        let (dir, nodes, names) = tax_fixture();
        let merged_path = dir.path().join("merged.dmp");
        let mut m = std::fs::File::create(&merged_path).unwrap();
        writeln!(m, "99999\t|\t28211\t|").unwrap();
        let tax = Taxonomy::from_dump(&nodes, &names)
            .unwrap()
            .with_merged(&merged_path)
            .unwrap();
        assert_eq!(tax.resolve(99999), 28211);
        assert_eq!(tax.resolve(28211), 28211);
    }

    #[test]
    fn unknown_tax_id_errors() {
        let (_d, nodes, names) = tax_fixture();
        let tax = Taxonomy::from_dump(&nodes, &names).unwrap();
        assert!(matches!(tax.lineage(99999), Err(TaxonomyError::Unknown(_))));
    }

    #[test]
    fn lineage_includes_self_first_root_last() {
        let (_d, nodes, names) = tax_fixture();
        let tax = Taxonomy::from_dump(&nodes, &names).unwrap();
        let p = tax.lineage(2).unwrap();
        assert_eq!(p, vec![2, 1]);
    }
}
