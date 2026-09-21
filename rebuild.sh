#!/usr/bin/env bash
# Regenerate tables, figures, and the compiled report from scratch.
# Usage: ./rebuild.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT/3_code"

echo "== 1/3 aggregating run results into tables =="
python scripts/aggregate.py

echo "== 2/3 generating figures =="
python scripts/make_figures.py

echo "== 3/3 compiling report =="
cd "$ROOT/5_report"

# TODO: no report compiler has been chosen yet. Pick one and replace this
# block, e.g.:
#   latexmk -pdf -cd template/main.tex -outdir=../5_report/build
#   latexmk -xelatex -cd template/main.tex -outdir=../5_report/build
#   pandoc chapters/*.md -o build/thesis.pdf --template=template/thesis.latex
echo "!! Report compile step not configured yet. Edit rebuild.sh once you pick a toolchain (latexmk/xelatex/pandoc)." >&2
exit 1
