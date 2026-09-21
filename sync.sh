#!/usr/bin/env bash
# Commit local changes, then sync report/ with Overleaf and the whole repo with GitHub.
# Usage: ./sync.sh ["commit message"]   (asks for a message if there are changes and none is given)
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

if [ -n "$(git status --porcelain)" ]; then
    git status --short
    msg="${1:-}"
    if [ -z "$msg" ]; then
        read -rp "Commit message: " msg
    fi
    if [ -z "$msg" ]; then
        echo "!! Empty commit message, nothing done." >&2
        exit 1
    fi
    echo "== 0/4 commit =="
    git add -A
    git commit -m "$msg"
fi

echo "== 1/4 GitHub -> local =="
git pull --no-rebase origin main

echo "== 2/4 Overleaf -> report/ =="
git subtree pull --prefix=report overleaf main --squash -m "Merge Overleaf changes into report/"

echo "== 3/4 report/ -> Overleaf =="
git subtree push --prefix=report overleaf main

echo "== 4/4 local -> GitHub =="
git push origin main
