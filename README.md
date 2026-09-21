# Look Where It Changed

Thesis: robot-guided active vision for updating Gaussian Splatting scenes.

## Folders

```
thesis/
├── admin/     forms, university guides, submissions/
├── papers/    related_work.xlsx; pdf/ is Zotero's linked-files folder (not in git)
├── code/      point_simulation/ — next-best-view simulator (see its notes.md)
├── data/      datasets (not in git)
├── runs/      experiment outputs (not in git)
├── report/    Overleaf project, synced as a git subtree
│   └── figures/auto/ written by code, figures/manual/ hand-made (keep sources)
├── archive/   dead work, local only (not in git)
└── sync.sh
```

## Sync

One GitHub repo (`origin`) holds everything tracked. Only `report/` is pushed to
Overleaf (`overleaf` remote) via `git subtree`.

```bash
./sync.sh "what changed"
```

`sync.sh` commits any local changes (asks for a message if none is given), pulls GitHub, pulls Overleaf into `report/`, pushes `report/` to
Overleaf, then pushes to GitHub. Run it before editing `.tex` files locally.

## Zotero

`papers/pdf/` is Zotero's linked attachment base directory: never move or
rename it without updating Zotero first. The bibliography goes to
`report/references.bib` via a Better BibTeX auto-export.
