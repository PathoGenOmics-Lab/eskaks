#!/usr/bin/env bash
# Repo hygiene check for ADDED lines only, so it never trips on pre-existing content.
# Fails if a change introduces:
#   - an em-dash (project preference is to avoid them; use a hyphen, comma, or parens)
#   - a leftover dbg!() macro
#   - a merge-conflict marker
#
# Usage:
#   scripts/check-hygiene.sh              # staged diff (used by the pre-commit hook)
#   scripts/check-hygiene.sh RANGE        # a git range, e.g. origin/main...HEAD (CI)
set -eu

range="${1:---cached}"

# Added lines of the diff (body lines starting with a single '+', never '+++').
# Excludes: generated golden fixtures (data, not prose), the example reports (also
# generated: the binary writes them, and the prose inside them comes from
# src/report/, which this check already covers at its source), and the hygiene
# tooling itself, whose source necessarily spells out the very patterns it forbids.
added="$(git diff --no-color -U0 $range -- . \
  ':(exclude)tests/golden/**' \
  ':(exclude)docs/assets/example-report*.html' \
  ':(exclude)scripts/check-hygiene.sh' \
  ':(exclude).githooks/**' 2>/dev/null \
  | grep -E '^\+' | grep -Ev '^\+\+\+' || true)"

# Build the em-dash byte sequence at runtime so this script contains none itself
# (and so we avoid grep -P, which BSD grep lacks).
emdash="$(printf '\xe2\x80\x94')"

fail=0
check() { # $1=description  $2=grep flag (-F|-E)  $3=pattern
  hits="$(printf '%s\n' "$added" | grep -n "$2" -- "$3" || true)"
  if [ -n "$hits" ]; then
    printf 'hygiene: %s\n' "$1" >&2
    printf '%s\n' "$hits" | sed 's/^/    /' >&2
    fail=1
  fi
}

check "an added line contains an em-dash (use a hyphen, comma, or parentheses)" -F "$emdash"
check "an added line contains a leftover dbg!() macro" -E 'dbg!\('
check "an added line contains a merge-conflict marker" -E '^\+(<<<<<<<|>>>>>>>)'

if [ "$fail" -ne 0 ]; then
  echo "hygiene check failed. Fix the above, or bypass with 'git commit --no-verify' if intentional." >&2
  exit 1
fi
echo "hygiene: OK"
