#!/usr/bin/env bash
# perfgate — the rsomics-* "must outperform upstream" release gate.
#
# Runs ours vs an upstream binary on one fixture under hyperfine, records
# full provenance (commit, upstream version, machine, fixture sha256,
# timing distribution) to .autopilot/state/perf-<date>.md, and exits
# non-zero unless ours is strictly faster. Equal-to-upstream is a failure
# by contract: a Rust port that only matches its C reference is ecosystem
# noise.
#
#   scripts/perfgate.sh \
#     --name fasta-stats-st \
#     --fixture     /path/big.fa \
#     --ours-bin    target/release/rsomics-fasta-stats --ours-args '-t 1 FIX' \
#     --upstream-bin seqkit --upstream-args 'stats -j 1 FIX' \
#     --upstream-version 'seqkit version'
#
# In --ours-args / --upstream-args the bare words FIX and OUT are replaced
# by the input fixture and a per-side scratch output path. bin + fixture
# are symlinked into a space-free tmp dir before timing, so substituted
# paths never contain spaces and this works from the apostrophe/space repo
# path under --shell=none. A run that does not clear >1.0x still gets
# written (FAIL verdict) — a recorded miss is the signal that an
# optimisation pass is owed before any tag.
set -euo pipefail

NAME=  FIXTURE=  FIXTURE2=  OURS_BIN=  OURS_ARGS=  UP_BIN=  UP_ARGS=  UP_VER=  WARMUP=3  MINRUNS=10
while [ $# -gt 0 ]; do
  case "$1" in
    --name)             NAME=$2; shift 2 ;;
    --fixture)          FIXTURE=$2; shift 2 ;;
    --fixture2)         FIXTURE2=$2; shift 2 ;;
    --ours-bin)         OURS_BIN=$2; shift 2 ;;
    --ours-args)        OURS_ARGS=$2; shift 2 ;;
    --upstream-bin)     UP_BIN=$2; shift 2 ;;
    --upstream-args)    UP_ARGS=$2; shift 2 ;;
    --upstream-version) UP_VER=$2; shift 2 ;;
    --warmup)           WARMUP=$2; shift 2 ;;
    --min-runs)         MINRUNS=$2; shift 2 ;;
    *) echo "perfgate: unknown arg $1" >&2; exit 2 ;;
  esac
done
[ -n "$NAME" ] && [ -n "$FIXTURE" ] && [ -n "$OURS_BIN" ] && [ -n "$UP_BIN" ] \
  || { echo "perfgate: --name --fixture --ours-bin --upstream-bin required" >&2; exit 2; }
[ -f "$FIXTURE" ] || { echo "perfgate: fixture not found: $FIXTURE" >&2; exit 2; }
command -v "$UP_BIN" >/dev/null || [ -x "$UP_BIN" ] \
  || { echo "perfgate: upstream not found: $UP_BIN" >&2; exit 2; }

repo_root=$(git rev-parse --show-toplevel); cd "$repo_root"

work=$(mktemp -d /tmp/perfgate.XXXXXX)
trap 'rm -rf "$work"' EXIT
ln -s "$(cd "$(dirname "$FIXTURE")" && pwd)/$(basename "$FIXTURE")" "$work/fixture"
if [ -n "$FIXTURE2" ]; then
  [ -f "$FIXTURE2" ] || { echo "perfgate: fixture2 not found: $FIXTURE2" >&2; exit 2; }
  ln -s "$(cd "$(dirname "$FIXTURE2")" && pwd)/$(basename "$FIXTURE2")" "$work/fixture2"
fi
ln -s "$(cd "$(dirname "$OURS_BIN")" && pwd)/$(basename "$OURS_BIN")" "$work/ours"
up_abs=$(command -v "$UP_BIN" || echo "$UP_BIN")
ln -s "$up_abs" "$work/upstream"

sha256() { if command -v sha256sum >/dev/null; then sha256sum "$1" | awk '{print $1}';
           else shasum -a 256 "$1" | awk '{print $1}'; fi; }
cpu_id() { if [ "$(uname)" = "Darwin" ]; then sysctl -n machdep.cpu.brand_string;
           else awk -F: '/model name/{print $2; exit}' /proc/cpuinfo | sed 's/^ //'; fi; }
ncores() { if [ "$(uname)" = "Darwin" ]; then sysctl -n hw.ncpu; else nproc; fi; }

commit=$(git rev-parse --short HEAD)
dirty=$(git diff --quiet || echo "+dirty")
fix_size=$(wc -c < "$FIXTURE" | tr -d ' ')
fix_sha=$(sha256 "$FIXTURE")
up_ver_str=$([ -n "$UP_VER" ] && eval "$UP_VER" 2>&1 | head -1 || echo "n/a")

# FIX2 before FIX (FIX is a prefix of FIX2).
subst() { local s=$1; s=${s//FIX2/$work/fixture2}; s=${s//FIX/$work/fixture}; s=${s//OUT/$work/$2}; echo "$s"; }
ours_args=$(subst "$OURS_ARGS" out.ours)
up_args=$(subst "$UP_ARGS" out.up)

tmp_json=$(mktemp); trap 'rm -rf "$work" "$tmp_json"' EXIT
hyperfine --warmup "$WARMUP" --min-runs "$MINRUNS" --shell=none \
  --export-json "$tmp_json" \
  -n ours     "$work/ours $ours_args" \
  -n upstream "$work/upstream $up_args" >/dev/null

read -r o_mean o_sd o_min o_max u_mean u_sd u_min u_max ratio verdict <<EOF
$(python3 - "$tmp_json" <<'PY'
import json, sys
r = json.load(open(sys.argv[1]))["results"]
o = next(x for x in r if x["command"] == "ours")
u = next(x for x in r if x["command"] == "upstream")
ratio = u["mean"] / o["mean"]
print(o["mean"], o["stddev"], o["min"], o["max"],
      u["mean"], u["stddev"], u["min"], u["max"],
      f"{ratio:.4f}", "PASS" if ratio > 1.0 else "FAIL")
PY
)
EOF

out=".autopilot/state/perf-$(date +%Y-%m-%d).md"
mkdir -p "$(dirname "$out")"
{
  echo "## $NAME — $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo
  echo "- ours: \`$(basename "$OURS_BIN") $OURS_ARGS <fixture>\` @ ${commit}${dirty}"
  echo "- upstream: \`$UP_BIN $UP_ARGS <fixture>\` — version: ${up_ver_str}"
  echo "- machine: $(uname -sm) | $(cpu_id) | $(ncores) cores"
  echo "- fixture: $FIXTURE — ${fix_size} bytes — sha256 ${fix_sha}"
  [ -n "$FIXTURE2" ] && echo "- fixture2: $FIXTURE2 — $(wc -c < "$FIXTURE2" | tr -d ' ') bytes — sha256 $(sha256 "$FIXTURE2")"
  echo "- hyperfine: warmup ${WARMUP}, min-runs ${MINRUNS}, shell=none"
  echo
  printf '| side | mean (s) | σ | min | max |\n|---|---|---|---|---|\n'
  printf '| ours | %.4f | %.4f | %.4f | %.4f |\n' "$o_mean" "$o_sd" "$o_min" "$o_max"
  printf '| upstream | %.4f | %.4f | %.4f | %.4f |\n' "$u_mean" "$u_sd" "$u_min" "$u_max"
  echo
  echo "**ratio (upstream/ours): ${ratio}× → ${verdict}** (contract: strictly > 1.0×)"
  echo
} >> "$out"

echo "perfgate $NAME: ${ratio}x → $verdict  (recorded → $out)"
[ "$verdict" = "PASS" ]
