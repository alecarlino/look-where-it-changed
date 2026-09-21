# Thesis

Computer vision thesis: code, data, experiment runs, and report live here in
one repo, wired so everything downstream of raw data can be regenerated with
one command (`./rebuild.sh`).

## Folder map

```
thesis/
├── 0_admin/            deadlines, forms, meeting notes, submissions/
├── 1_papers/            pdf/, references.bib, notes/
├── 2_data/              README.md, raw/, processed/, splits/
├── 3_code/
│   ├── paths.py          every path defined once (ROOT, DATA, RUNS, FIGURES, TABLES)
│   ├── viz.py            plot style + save_fig()
│   ├── tables.py         save_table() -> .md and .tex
│   ├── src/              my code: data/ models/ train/ eval/
│   ├── vendor/           third-party code, one folder each, never edited (see ORIGIN.md in each)
│   ├── scripts/          prepare_data.py, train.py, evaluate.py, aggregate.py, make_figures.py
│   ├── configs/          baseline.yaml, ablations/
│   ├── nbv_simulation/   next best view simulator over point clouds
│   │                     (see its notes.md for the technical part)
│   ├── sandbox/          scratch, dated folders, never imported by src/
│   └── tests/
├── 4_runs/              baseline/, ablations/<axis>/<variant>/, each with
│                        config.yaml, metrics.csv, log.txt, ckpt/
├── 5_report/
│   ├── chapters/         .md and .tex chapters coexist
│   ├── figures/auto/     written by code only (via viz.save_fig())
│   ├── figures/manual/   hand-made figures, keep .svg sources
│   ├── tables/           generated, .md + .tex per table (via tables.save_table())
│   ├── template/
│   └── build/            compiled output (gitignored)
├── 9_archive/           dead work, dated folders, _unsorted/ for anything unclassified
├── rebuild.sh           figures + tables + report pdf in one command
└── .gitignore
```

## Data flow (one-way)

```
2_data  →  3_code  →  4_runs  →  5_report
(raw)      (train/      (results)   (figures, tables,
 eval)                    chapters, pdf)
```

Nothing upstream ever reads from something downstream: `3_code` never writes
into `2_data/raw`, `4_runs` never edits `3_code`, and `5_report` only ever
*reads* from `4_runs` via `3_code/scripts/aggregate.py` and
`3_code/scripts/make_figures.py`. Every path used anywhere is defined once in
[`3_code/paths.py`](3_code/paths.py) — no hardcoded absolute paths elsewhere.

## Regenerating everything

```bash
./rebuild.sh
```

Runs, in order: `aggregate.py` (collect `4_runs/*/metrics.csv` into
`5_report/tables/`) → `make_figures.py` (regenerate `5_report/figures/auto/`)
→ compile the report. The compile step isn't wired up yet — see the TODO in
`rebuild.sh`.

## Status

This structure was scaffolded fresh (the folder was empty when set up on
2026-09-14) — `src/`, `scripts/`, `vendor/`, data, and runs are all still to
be filled in.
