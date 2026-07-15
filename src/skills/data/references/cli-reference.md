# CLI reference — `sqanti_browser`

Single flat command (no subcommands). Source of truth: the argparse block in
`sqanti_browser.py`. Run as `python -m sqanti_browser …` (recommended), `sqanti_browser …`
(after install), or `python sqanti_browser.py …`.

## Required

| Flag | Value | Notes |
|------|-------|-------|
| `--gtf` | path | SQANTI3 corrected GTF (`*_corrected.gtf`). Missing file → exit 1. |
| `--classification` | path | SQANTI3 classification table (`*_classification.txt`). Missing → exit 1. |
| `--output` | dir | Output directory (created with parents). Becomes the hub root. |
| `--genome` | name | Assembly name (`hg38`, `mm10`, or any custom name with `--twobit`). Used as the genome subdir and track-file prefix. |

## Optional inputs / tracks

| Flag (dest) | Value | Effect |
|------|-------|--------|
| `--chrom-sizes` | path | Explicit chrom.sizes; otherwise derived. |
| `--twobit` | `.2bit` | Non-reference genome: computes chrom.sizes via `twoBitInfo`, copies the 2bit into the hub, writes a genome HTML + `defaultPos`. Requires `twoBitInfo`. |
| `--star-sj` | `SJ.out.tab` | STAR splice-junction track → `<genome>_star_sj.bb` (strand 1→`+`, 2→`-`). |
| `--CAGE-peak` (`cage_peak`) | BED | TSS validation track → `<genome>_cage_peaks.bb`. |
| `--polyA-peak` (`polya_peak`) | BED | TTS validation track → `<genome>_polya_peaks.bb`. |
| `--refGTF` (`ref_gtf`) | GTF | Reference annotation for direct comparison → `<genome>_reference.bb` (gray). Recommended. |

Peak/reference/chrom inputs are filtered to chromosomes and coordinate bounds present in
the genome; out-of-range features are dropped (with warnings).

## Display / behavior

| Flag | Default | Effect |
|------|---------|--------|
| `--tables` | off | After the hub, also generate interactive HTML tables per category into `<output>/table_reports/`. |
| `--sort-by CHOICE` | `basic` | Isoform order within genomic position. Choices below. |
| `--no-category-tracks` | off | Main track only; skip per-category tracks (faster). |
| `--category-tracks LIST` | all | Comma-separated abbreviations to limit per-category tracks. Unknown abbrevs warn and are skipped. |
| `--no-highlight` | off | Disable highlight color for the top-FL isoform per group. |
| `--hub-name NAME` | (hub) | Display name + prefix on every track label. Use distinct names to compare multiple hubs. |
| `--my-palette JSON` | defaults | Custom colors. Missing file → exit 1. See custom-coloring.md. |
| `--trix` / `--no-trix` | trix on | `--no-trix` skips `trix.ix`/`trix.ixx`. `--trix` is the default; needs `ixIxx` (silently skipped if absent). Mutually exclusive. |

### `--sort-by` choices

`basic` (genomic position; pipeline order for ties — default), `none` (genomic position;
GTF-file order for ties), `iso_exp`, `length`, `FL`, `diff_to_TSS`, `diff_to_TTS`,
`diff_to_gene_TSS`, `diff_to_gene_TTS`, `dist_to_CAGE_peak`, `dist_to_polyA_site`.
Distance metrics sort ascending (closest first); the rest descending (largest first).
(Note: the wiki's CLI table mislabels the default as `none`; the code default is `basic`.)

### `--category-tracks` abbreviations

`FSM`→full-splice_match, `ISM`→incomplete-splice_match, `NIC`→novel_in_catalog,
`NNC`→novel_not_in_catalog, `antisense`, `genic_intron`, `genic_genomic`(or `genic`)→genic,
`intergenic`, `fusion`. Example: `--category-tracks FSM,ISM,NIC`.

## Debug / info

| Flag | Effect |
|------|--------|
| `--validate-only` | Validate tools + inputs, then exit (fast environment check). |
| `--dry-run` | Build the intermediate BED (classification merged) and exit before bigBed/hub. |
| `--keep-temp` | Keep the temp dir (also kept automatically on failure). |
| `--version` | Print `sqanti_browser 1.1.2` and exit. |

> **`--install-only` is not a `sqanti_browser` flag** — it belongs to the test harness
> (`python tests/test_sqanti_browser.py --install-only`). For a CLI env check use `--validate-only`.

## Standalone report CLI — `src/filter_isoforms.py`

Generate the interactive HTML tables **without** building a hub:

```bash
python src/filter_isoforms.py --classification classification.txt --output-dir reports [--include-sequences]
```

- `--classification` (required), `--output-dir` (default `isoform_reports`),
  `--include-sequences` (keep the `ORF_seq` column; off by default).
