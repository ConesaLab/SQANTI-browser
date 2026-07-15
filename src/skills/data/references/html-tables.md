# Interactive HTML tables

Add `--tables` to a run to emit filterable HTML reports of the classification data —
explorable in any browser, no UCSC hub or account required.

```bash
python -m sqanti_browser \
    --gtf corrected.gtf --classification classification.txt \
    --output my_hub --genome hg38 --tables
```

## Output

Reports land in `<output>/table_reports/`:

- `complete_transcriptome_isoforms.html` — every isoform, all categories
- `<category>_isoforms.html` — one per structural category, e.g.
  `full-splice_match_isoforms.html`, `incomplete-splice_match_isoforms.html`,
  `novel_in_catalog_isoforms.html`, `novel_not_in_catalog_isoforms.html`,
  `genic_isoforms.html`, `antisense_isoforms.html`, `fusion_isoforms.html`,
  `intergenic_isoforms.html`, `genic_intron_isoforms.html`

Each file carries all ~48 classification metrics for its isoforms. The `ORF_seq` column is
excluded by default (huge); re-add it with `--include-sequences` on the standalone CLI.

## Internet / CDN requirement

Tables are built on DataTables/jQuery loaded from a **CDN**, so the page needs internet
access when opened in the browser to render fully. (The files themselves are portable, but
rendering pulls the JS from the CDN.)

## In-page features

- **Column filters**: numeric ranges (`100:1000`, `100:`, `:1000`; a bare number = exact),
  dropdowns on categorical fields, substring text search (case-insensitive). A global search
  box filters across all columns.
- **Sortable columns**: click a header (first click ascending, second descending).
- **Row selection**: click rows to toggle multi-select (highlighted).
- **Export**: download the current (filtered/selected) rows to CSV.
- **Generate Filter String**: produces a newline-separated list of isoform IDs (the
  `isoform` column) — selected rows if any, otherwise the full filtered set. Use these IDs
  in the UCSC Table Browser to build a custom subset track (see
  [subset-sessions.md](subset-sessions.md)).
- **Generate Trix String**: builds a UCSC search string. Filter-based (from current
  dropdown filters) or row-based (from a selected isoform). Paste into the UCSC search box to
  jump to matching isoforms. See [filtering-and-search.md](filtering-and-search.md) for token
  syntax.

## Standalone generation (no hub)

Generate the same reports without building a hub:

```bash
python src/filter_isoforms.py \
    --classification classification.txt \
    --output-dir reports
```

- `--classification` (required), `--output-dir` (default `isoform_reports`),
  `--include-sequences` (keep the `ORF_seq` column, makes files much larger).

## Limitations

- No genome context — use the UCSC hub for that.
- Trix strings cover dropdown/categorical filters only; numeric ranges are handled by UCSC
  track filters, not Trix.
- Very large datasets (100K+ isoforms) load slowly; filter or split by category first.
