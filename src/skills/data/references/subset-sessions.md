# Curated subset sessions

Build a UCSC custom track containing only a hand-picked set of isoforms, without rebuilding
the hub. The flow: filter/select in the HTML report → export isoform IDs → paste into the
UCSC Table Browser → output a custom track. The full hub stays available for context.

## Workflow

| Step | Where | Action |
|------|-------|--------|
| 1 | SQANTI-browser | Run with `--tables` to get HTML reports |
| 2 | HTML report | Filter columns and/or select rows |
| 3 | HTML report | Click **Generate Filter String** → copy isoform IDs |
| 4 | UCSC Table Browser | Paste IDs into Identifiers, output = custom track |
| 5 | UCSC Genome Browser | View the subset alongside the hub |

## Step 1 — generate tables

```bash
python -m sqanti_browser \
    --gtf your_corrected.gtf --classification your_classification.txt \
    --output my_hub --genome hg38 --tables
```

Reports go to `<output>/table_reports/`: `complete_transcriptome_isoforms.html` (all
categories) and `<category>_isoforms.html` per category. You can also generate tables alone:
`python src/filter_isoforms.py --classification your_classification.txt --output-dir table_reports`.

## Step 2 — curate in the report

Filter with column headers: numeric ranges (`100:1000`, `100:`, `:1000`, bare number =
exact), dropdown exact-match, substring text search. Category reports are pre-restricted to
one structural category; use the complete report to work across categories. Optionally click
rows to toggle multi-select (highlighted).

## Step 3 — export IDs

Click **Generate Filter String**. Behavior:

- One or more rows selected → exports **only selected** isoforms.
- No rows selected → exports **all rows matching current filters**.

The dialog lists isoform IDs, one per line (from the `isoform` column). Copy them all.

**Critical:** these IDs must match the bigBed track's `name` field / the `isoform` values in
`classification.txt` **exactly** — paste unchanged, no edits.

## Step 4 — Table Browser custom track

With the hub connected in UCSC:

1. **Tools → Table Browser**.
2. Select the same **genome assembly** as the hub (e.g. Human / hg38).
3. **group** → your Track Hub group.
4. **track** → the SQANTI3 bigBed track that contains those IDs (main track, or the matching
   category track).
5. **region** → genome (or a locus).
6. Click **Identifiers**, paste the ID list (one per line) for the `name` field.
7. **output format** → **custom track** (or GTF/BED to download).
8. Run → **get custom track in Genome Browser**.

UCSC returns a temporary custom track of only those isoforms, on top of the full hub.

## Troubleshooting

| Issue | Fix |
|-------|-----|
| No IDs in dialog | Apply filters or select rows before Generate Filter String. |
| Table Browser finds no matches | Hub must be connected; IDs must exactly equal `isoform` values in `classification.txt`. |
| Subset track empty | Pick the bigBed track that actually contains those IDs. |
| Very large ID list | UCSC limits custom track size — split into batches or narrow filters first. |
| Want search instead of a fixed list | Use Generate Trix String + hub search ([filtering-and-search.md](filtering-and-search.md)). |
