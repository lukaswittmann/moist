#!/bin/bash
# sync-version.sh -- propagate VERSION to all version declarations
#
# VERSION may carry an optional SemVer pre-release suffix, e.g.
#   0.6.0         (stable)
#   0.6.0-alpha.1 (alpha)
#
# Usage: config/sync-version.sh [--check]
#   --check   exit 1 if any file would change, without writing (for CI)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VERSION_FILE="$REPO_ROOT/VERSION"

if [ ! -f "$VERSION_FILE" ]; then
   echo "error: VERSION file not found at $VERSION_FILE" >&2
   exit 1
fi

FULL="$(tr -d '[:space:]' < "$VERSION_FILE")"

semver='^[0-9]+\.[0-9]+\.[0-9]+(-[0-9A-Za-z.]+)?$'
if ! [[ $FULL =~ $semver ]]; then
   echo "error: VERSION '$FULL' is not MAJOR.MINOR.PATCH[-PRERELEASE]" >&2
   exit 1
fi

# Numeric base (strip any -PRERELEASE) for resolver-safe targets.
BASE="${FULL%%-*}"
IFS='.' read -r MAJOR MINOR PATCH <<< "$BASE"

# PEP 440 spelling of the pre-release: alpha.1 -> a1, beta -> b, rc.2 -> rc2
PRE="${FULL#"$BASE"}"
PRE="${PRE#-}"
num="${PRE##*.}"
if [ "$num" = "$PRE" ]; then num=""; fi   # no numeric segment after a dot
case "${PRE%%.*}" in
   alpha|a) PEP440="${BASE}a${num}" ;;
   beta|b)  PEP440="${BASE}b${num}" ;;
   rc|c)    PEP440="${BASE}rc${num}" ;;
   *)       PEP440="$BASE" ;;
esac

CHECK=0
if [ "${1:-}" = "--check" ]; then
   CHECK=1
fi
CHANGED_FILES=()

tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

# Apply sed expressions to a file; with --check only record that it would
# change. Writing through a redirect keeps the file's mode.
sub() {
   local file=$1
   shift
   sed "$@" "$file" > "$tmp"
   if cmp -s "$tmp" "$file"; then
      return 0
   fi
   if [ "$CHECK" -eq 1 ]; then
      CHANGED_FILES+=("$file")
   else
      cat "$tmp" > "$file"
   fi
}

echo "Syncing version $FULL (base $BASE, wheel $PEP440) to all targets..."

# meson.build: the full string (meson tolerates the suffix); python/meson.build:
# the PEP 440 one. Only the project line: its version starts with a digit,
# dependency constraints start with >=.
project_version="s/^\(  version: '\)[0-9][^']*'"
sub "$REPO_ROOT/meson.build" -e "$project_version/\1$FULL'/"
sub "$REPO_ROOT/python/meson.build" -e "$project_version/\1$PEP440'/"

# python/pyproject.toml declares the version dynamic (meson-python reads
# python/meson.build) and moist/__init__.py asks the library: nothing to sync.

# version.f90: display string (full) and compact array (base integers)
sub "$REPO_ROOT/src/moist/version.f90" \
   -e "s/moist_version_string = \"[^\"]*\"/moist_version_string = \"$FULL\"/" \
   -e "s/moist_version_compact(3) = \[[0-9, ]*\]/moist_version_compact(3) = [$MAJOR, $MINOR, $PATCH]/"

if [ "$CHECK" -eq 0 ]; then
   echo "Done."
elif [ "${#CHANGED_FILES[@]}" -gt 0 ]; then
   echo "" >&2
   echo "error: version files are out of sync with VERSION." >&2
   echo "Run 'config/sync-version.sh' and commit the result." >&2
   printf '  %s\n' "${CHANGED_FILES[@]}" >&2
   exit 1
else
   echo "All version files are in sync."
fi
