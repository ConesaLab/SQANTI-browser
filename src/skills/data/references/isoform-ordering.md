# Isoform ordering

`--sort-by` controls the order isoforms appear **within each genomic position**. It only
affects display order in the browser — it never filters or removes isoforms.

```bash
python -m sqanti_browser ... --sort-by iso_exp
```

## Why ordering is only local

bigBed (required by UCSC) must be sorted by genomic position first. You cannot globally
reorder the whole track by a metric — a highly-expressed isoform that starts further
downstream still appears after everything upstream of it. `--sort-by` decides ordering
inside a position group only.

## Three-level hierarchy

1. **Genomic position** — chromosome then start coordinate. Fixed by bigBed, cannot change.
2. **Reference transcript group** — isoforms sharing a reference are grouped together
   (automatic).
3. **Tie-breaker within a group** — set by `--sort-by`.

## `--sort-by` choices

Default is **`basic`** (not `none`).

| Value | Within-group order |
|-------|--------------------|
| `basic` (default) | Pipeline order |
| `none` | Preserves GTF file order — matches what you'd see loading the GTF directly into UCSC |
| `iso_exp` | Isoform expression, highest first |
| `FL` | Full-length read count, highest first |
| `length` | Transcript length, longest first |
| `diff_to_TSS` | Distance to ref TSS, closest first |
| `diff_to_TTS` | Distance to ref TTS, closest first |
| `diff_to_gene_TSS` / `diff_to_gene_TTS` | Distance to gene TSS/TTS, closest first |
| `dist_to_CAGE_peak` | Distance to CAGE peak, closest first |
| `dist_to_polyA_site` | Distance to polyA site, closest first |

**Direction rule:** distance metrics (`diff_to_*`, `dist_to_*`) sort **ascending** (closest
first); all other metrics sort **descending** (largest/highest first).

## NA handling

If the chosen sort column has NA values: a warning is printed, NA rows sort last, and if the
whole column is NA or sorting fails the tool falls back to `none` (GTF order).

## Examples

```bash
# Highest-expressed isoform first within each position
python -m sqanti_browser ... --sort-by iso_exp

# Prioritize CAGE-validated TSS (requires --CAGE-peak input)
python -m sqanti_browser ... --CAGE-peak CAGE_peaks.bed --sort-by dist_to_CAGE_peak

# Reproduce a direct GTF upload's order
python -m sqanti_browser ... --sort-by none
```

Note: the wiki's command_line_reference table lists the default as `none` — that is wrong;
the code default is `basic`. See [cli-reference.md](cli-reference.md).
