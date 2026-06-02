#!/bin/bash
# sync-version.sh -- propagate VERSION to all version declarations
#
# VERSION may carry an optional SemVer pre-release suffix, e.g.
#   0.6.0         (stable)
#   0.6.0-alpha.1 (alpha)
#
# Usage: config/sync-version.sh [--check]
#   --check   exit 1 if any file would change (for CI validation)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VERSION_FILE="$REPO_ROOT/VERSION"

if [ ! -f "$VERSION_FILE" ]; then
   echo "error: VERSION file not found at $VERSION_FILE" >&2
   exit 1
fi

FULL="$(tr -d '[:space:]' < "$VERSION_FILE")"

# Validate SemVer with optional pre-release: MAJOR.MINOR.PATCH[-PRERELEASE]
if ! echo "$FULL" | grep -qE '^[0-9]+\.[0-9]+\.[0-9]+(-[0-9A-Za-z.]+)?$'; then
   echo "error: VERSION '$FULL' is not MAJOR.MINOR.PATCH[-PRERELEASE]" >&2
   exit 1
fi

# Numeric base (strip any -PRERELEASE) for resolver-safe targets.
BASE="${FULL%%-*}"
IFS='.' read -r MAJOR MINOR PATCH <<< "$BASE"

# Translate a SemVer pre-release (e.g. "alpha.1") to PEP 440 spelling
PRE="${FULL#"$BASE"}"
PEP440="$BASE"
if [ -n "$PRE" ]; then
   label="${PRE#-}"          # alpha.1
   name="${label%%.*}"       # alpha
   num="${label##*.}"        # 1
   if [ "$num" = "$label" ]; then num=""; fi   # no numeric segment after the dot
   case "$name" in
      alpha|a) PEP440="${BASE}a${num}" ;;
      beta|b)  PEP440="${BASE}b${num}" ;;
      rc|c)    PEP440="${BASE}rc${num}" ;;
      *)       PEP440="$BASE" ;;
   esac
fi

# Portable in-place edit
sub() {
   local expr=$1 file=$2 tmp
   tmp="$(mktemp "${file}.XXXXXX")"
   if sed "$expr" "$file" > "$tmp"; then
      mv "$tmp" "$file"
   else
      rm -f "$tmp"
      return 1
   fi
}

# Files that will be checked in --check mode
TARGETS=(
   meson.build
   fpm.toml
   python/meson.build
   python/pyproject.toml
   src/moist/version.f90
   doc/src/version.f90
)

echo "Syncing version $FULL (base $BASE, wheel $PEP440) to all targets..."

# 1. meson.build (root) -- version: 'X.Y.Z[-pre]' (full string; meson tolerates it)
#    Match only lines where version value starts with a digit (skip >=... patterns)
sub "s/^\(  version: '\)[0-9][^']*'/\1$FULL'/" "$REPO_ROOT/meson.build"

# 2. fpm.toml -- version = "X.Y.Z" (numeric base only; fpm rejects suffixes)
sub "s/^version = \"[^\"]*\"/version = \"$BASE\"/" "$REPO_ROOT/fpm.toml"

# 3. python/meson.build -- version: 'X.Y.Z[aN]' (PEP 440; project line only, not dependency)
sub "s/^\(  version: '\)[0-9][^']*'/\1$PEP440'/" "$REPO_ROOT/python/meson.build"

# 4. python/pyproject.toml -- version = "X.Y.Z[aN]" (PEP 440 wheel version)
sub "s/^version = \"[^\"]*\"/version = \"$PEP440\"/" "$REPO_ROOT/python/pyproject.toml"

# 5 & 6. version.f90 (src + doc) -- display string (full) + compact array (base ints)
for f90 in "$REPO_ROOT/src/moist/version.f90" "$REPO_ROOT/doc/src/version.f90"; do
   if [ -f "$f90" ]; then
      sub "s/moist_version_string = \"[^\"]*\"/moist_version_string = \"$FULL\"/" "$f90"
      sub "s/moist_version_compact(3) = \[[0-9, ]*\]/moist_version_compact(3) = [$MAJOR, $MINOR, $PATCH]/" "$f90"
   fi
done

echo "Done."

# --check mode: fail if the sync produced changes in version-related files only
if [ "${1:-}" = "--check" ]; then
   CHANGED=0
   for t in "${TARGETS[@]}"; do
      if ! git -C "$REPO_ROOT" diff --quiet -- "$t" 2>/dev/null; then
         CHANGED=1
      fi
   done
   if [ "$CHANGED" -ne 0 ]; then
      echo "" >&2
      echo "error: version files are out of sync with VERSION." >&2
      echo "Run 'config/sync-version.sh' and commit the result." >&2
      git -C "$REPO_ROOT" diff --stat -- "${TARGETS[@]}" >&2
      exit 1
   fi
   echo "All version files are in sync."
fi
