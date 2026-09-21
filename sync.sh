#!/usr/bin/env bash
# Sync report/ with Overleaf and the whole repo with GitHub.
# Usage: ./sync.sh   (commit your local changes first)
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "!! Uncommitted changes: commit them first (git add -A && git commit -m ...)." >&2
    exit 1
fi

echo "== 1/4 GitHub -> local =="
git pull --no-rebase origin main

echo "== 2/4 Overleaf -> report/ =="
git subtree pull --prefix=report overleaf main --squash -m "Merge Overleaf changes into report/"

echo "== 3/4 report/ -> Overleaf =="
git subtree push --prefix=report overleaf main

echo "== 4/4 local -> GitHub =="
git push origin main
