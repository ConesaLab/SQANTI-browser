# Filtering and search in UCSC

Two complementary ways to narrow transcripts once the hub is loaded: hub-native track
filters (dropdowns, text, ranges) and Trix text search in the search box.

## Hub-native filters

Right-click the track → **Configure** (not plain track settings) → set filters → **Submit**.
Filter definitions come from `src/constants.py` (`ORDERED_FILTERS`), and which fields render
as dropdowns vs text comes from `DROPDOWN_COLS`.

### Dropdown filters (categorical)

Exact-match pick lists on: `structural_category`, `subcategory`, `coding`, `FSM_class`,
`all_canonical`, `RTS_stage`, `bite`, `predicted_NMD`, `within_CAGE_peak`,
`polyA_motif_found`.

### Text filters (wildcard)

Free text with `*` wildcard on fields like `strand`, `associated_gene`,
`associated_transcript`, `min_cov_pos` (e.g. `ENSG*`, `BRCA1`).

### Range filters (numeric, with default caps)

Sliders default to the min/max baked into `ORDERED_FILTERS`. Values beyond these caps
aren't reachable by the slider by default — widen the box if you need to. Notable defaults:

| Field | Default range |
|-------|---------------|
| `length` | 0 – 150,000 (150 kb) |
| `exons` | 0 – 400 |
| `iso_exp`, `gene_exp`, `FL` | 0 – 100,000 |
| `diff_to_TSS` / `diff_to_TTS` / `diff_to_gene_TSS` / `diff_to_gene_TTS` | -100,000 – 100,000 |
| `dist_to_CAGE_peak` / `dist_to_polyA_site` | -10,000 – 10,000 |
| `ORF_length` / `CDS_length` | 0 – 50,000 |
| `perc_A_downstream_TTS`, `ratio_exp` (0–1), `ratio_TSS` (0–10) | percent / ratio ranges |

If filters or values seem missing, add `udcTimeout=5` to the URL and reconnect the hub.

## Trix text search

Type tokens into the UCSC search box. Basic search matches an isoform ID (`PB.1234.1`),
gene name (`BRCA1`), or gene ID (`ENSG00000012048`).

### Token syntax — all underscores

Attribute tokens join the field and value with underscores; hyphens and spaces in values
are normalized to underscores. Examples:

```
structural_category_full_splice_match      # full-splice_match category
structural_category_incomplete_splice_match
structural_category_novel_in_catalog
structural_category_novel_not_in_catalog
strand_plus            strand_minus
coding_coding          coding_non_coding
subcategory_mono_exon
predicted_NMD_True     within_CAGE_peak_True
polyA_motif_found_True RTS_stage_True   bite_True
```

Combine tokens with spaces for **AND** logic (all must match):

```
BRCA1 structural_category_novel_in_catalog coding_coding
```

### Rules and limits

- Case-insensitive; partial matches work for gene names/IDs.
- AND only — no OR, no numeric ranges (use track range filters for numbers).
- The HTML tables' **Generate Trix String** button builds these tokens for you from the
  current filters or a selected row — copy and paste into the search box.
- If search returns nothing: check `trix.ix`/`trix.ixx` exist in the genome dir, verify
  `trackDb.txt` has `searchTrix`, and clear cache with `udcTimeout=5`.

See [html-tables.md](html-tables.md) and [cli-reference.md](cli-reference.md).
