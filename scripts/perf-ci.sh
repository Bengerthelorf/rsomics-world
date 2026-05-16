#!/usr/bin/env bash
# perf-ci — the publish-time enforcement of the >1.0x contract.
#
#   scripts/perf-ci.sh <crate>
#
# Foundation libraries (crates/foundation/*) have no binary/upstream →
# pass (notice). Layer-B tool crates are resolved against
# scripts/perf-manifest.toml: `exempt` passes with its reason; `spec`
# builds the release binary, generates the fixture(s), runs perfgate and
# propagates its exit code (non-zero unless strictly >1.0x); anything
# else (needs_spec / absent) BLOCKS — a port ships only once proven
# faster than its reference.
set -euo pipefail

CRATE=${1:?usage: perf-ci.sh <crate>}
root=$(git rev-parse --show-toplevel); cd "$root"

manifest=$(find crates -mindepth 2 -maxdepth 5 -path "*/$CRATE/Cargo.toml" | head -n 1)
[ -n "$manifest" ] || { echo "::error::perf-ci: no crate dir for $CRATE"; exit 1; }
case "$manifest" in
  *crates/foundation/*)
    echo "::notice::perf-ci: $CRATE is a foundation library — no binary/upstream, perf-exempt by nature."
    exit 0 ;;
esac

cls=$(python3 - "$CRATE" <<'PY'
import sys, pathlib
try:
    import tomllib            # py3.11+
except ModuleNotFoundError:
    import tomli as tomllib   # py3.10 (4090): pip install --user tomli
crate = sys.argv[1]
m = tomllib.loads(pathlib.Path("scripts/perf-manifest.toml").read_text())
if crate in m.get("exempt", {}):
    print("EXEMPT", m["exempt"][crate].get("reason", "")); raise SystemExit
if crate in m.get("spec", {}):
    s = m["spec"][crate]
    print("SPEC")
    for k in ("fixture", "fixture2", "ours_args", "upstream_bin",
              "upstream_args", "upstream_version"):
        if k in s:
            print(f"{k}\t{s[k]}")
    raise SystemExit
print("BLOCK")
PY
)
kind=$(printf '%s\n' "$cls" | head -n1 | cut -d' ' -f1)

case "$kind" in
  EXEMPT)
    echo "::notice::perf-ci: $CRATE perf-exempt — $(printf '%s' "$cls" | head -n1 | cut -d' ' -f2-)"
    exit 0 ;;
  BLOCK|"")
    echo "::error::perf-ci: $CRATE has no verified perfgate spec (needs_spec/absent). A Layer-B port may not publish until it is proven strictly >1.0x vs its upstream. Add a [spec.$CRATE] to scripts/perf-manifest.toml and let it pass."
    exit 1 ;;
esac

# Fixed key set → plain vars (macOS ships bash 3.2; no `declare -A`).
f_fixture=  f_fixture2=  f_ours_args=  f_upstream_bin=  f_upstream_args=  f_upstream_version=
while IFS=$'\t' read -r k v; do
  case "$k" in
    fixture)          f_fixture=$v ;;
    fixture2)         f_fixture2=$v ;;
    ours_args)        f_ours_args=$v ;;
    upstream_bin)     f_upstream_bin=$v ;;
    upstream_args)    f_upstream_args=$v ;;
    upstream_version) f_upstream_version=$v ;;
  esac
done <<< "$cls"

mkfix() { # "<kind>:<a>:<b>" dst [seed] -> path
  local spec=$1 dst=$2 seed=${3:-} IFS=:
  set -- $spec
  # A gz fixture must keep a .gz name so perfgate's extension-preserving
  # symlink lets extension-sniffing upstreams (fastp) detect compression.
  [ "$1" = fastqgz ] && dst="$dst.fastq.gz"
  python3 scripts/mkfixture.py "$1" "$dst" "$2" "$3" $seed >/dev/null
  echo "$dst"
}

work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
fix=$(mkfix "$f_fixture" "$work/fix")
tdir=${CARGO_TARGET_DIR:-target}
args=(--name "ci-${CRATE}" --fixture "$fix"
      --ours-bin "$tdir/release/${CRATE}" --ours-args "$f_ours_args"
      --upstream-bin "$f_upstream_bin" --upstream-args "$f_upstream_args")
# Distinct seed so a (FIX) and b (FIX2) differ for two-input tools.
[ -n "$f_fixture2" ] && args+=(--fixture2 "$(mkfix "$f_fixture2" "$work/fix2" 424242)")
[ -n "$f_upstream_version" ] && args+=(--upstream-version "$f_upstream_version")

cargo build --release -p "$CRATE"
exec scripts/perfgate.sh "${args[@]}"
