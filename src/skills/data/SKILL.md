---
name: sqanti-browser
description: Use whenever the user is working with SQANTI-browser — turning SQANTI3 QC output (a corrected GTF + classification file) into a UCSC Genome Browser track hub for visualizing and curating long-read transcriptomes. Covers the `sqanti_browser` command and all its flags, required/optional inputs, the track-hub output (hub.txt, genomes.txt, trackDb, bigBed, groups, trix search index, per-category and validation tracks), interactive HTML tables (`--tables`), hosting the hub and validating with hubCheck, UCSC filtering and Trix search, isoform ordering (`--sort-by`), custom color palettes (`--my-palette`), non-reference genomes (`--twobit`), curated subset sessions via the Table Browser, and SQANTI-reads multi-sample workflows. Trigger for any mention of SQANTI-browser, "SQANTI3 to UCSC", building a track hub / trackDb / bigBed from SQANTI3, or visualizing classified isoforms in the UCSC browser.
---

# SQANTI-browser (v1.1.2)

SQANTI-browser converts **SQANTI3 QC output** into a **UCSC Genome Browser track hub** so
long-read transcriptomes can be explored, filtered, searched, and curated interactively.
It is a **command-line tool** — a single flat command, no subcommands and no meaningful
Python API. Drive it by constructing the right `sqanti_browser` invocation.

This is a **router skill**: use the decision table below and read the relevant
`references/` file before writing commands, rather than guessing flag names.

## Running it

Three equivalent forms (from the repo root unless installed):

- `python -m sqanti_browser …` — **recommended** (uses the active/conda Python, avoids
  wrong-Python and architecture mismatches).
- `sqanti_browser …` — after `pip install -e .` / `pip install sqanti-browser`.
- `python sqanti_browser.py …` — needs the project root importable.

## Canonical example (verified, bundled data)

```bash
python -m sqanti_browser \
  --gtf example/SQANTI3_QC_output/example_corrected.gtf \
  --classification example/SQANTI3_QC_output/example_classification.txt \
  --output example_output \
  --genome hg38
```

Produces `example_output/hub.txt`, `genomes.txt`, and `hg38/…_sqanti3.bb` + trackDb etc.
Then host the output dir and load the **raw** `hub.txt` URL in UCSC (My Data → Track Hubs).

## Required arguments (always)

| Flag | Meaning |
|------|---------|
| `--gtf` | SQANTI3 corrected GTF (`*_corrected.gtf`). |
| `--classification` | SQANTI3 classification table (`*_classification.txt`). |
| `--output` | Output directory (created if missing) — becomes the hub root. |
| `--genome` | Assembly name (e.g. `hg38`, `mm10`); names the genome subdir + track prefixes. |

## Decision table — read the reference before acting

| If the user wants to… | Read |
|---|---|
| Know exactly what a flag does / see all flags, defaults, choices | [`references/cli-reference.md`](references/cli-reference.md) |
| Understand required/optional inputs and the hub files produced | [`references/inputs-outputs.md`](references/inputs-outputs.md) |
| Install the tool + UCSC command-line utilities | [`references/installation.md`](references/installation.md) |
| Put the hub online and load it in UCSC (raw URL, hubCheck) | [`references/hosting.md`](references/hosting.md) |
| Generate/interpret the interactive HTML tables (`--tables`) | [`references/html-tables.md`](references/html-tables.md) |
| Filter in UCSC or build Trix search terms | [`references/filtering-and-search.md`](references/filtering-and-search.md) |
| Control isoform order (`--sort-by`) | [`references/isoform-ordering.md`](references/isoform-ordering.md) |
| Customize colors (`--my-palette`, highlight) | [`references/custom-coloring.md`](references/custom-coloring.md) |
| Use a non-reference genome (`--twobit`) | [`references/non-reference-genomes.md`](references/non-reference-genomes.md) |
| Show only a curated subset of isoforms in UCSC | [`references/subset-sessions.md`](references/subset-sessions.md) |
| Visualize a multi-sample SQANTI-reads experiment | [`references/sqanti-reads-integration.md`](references/sqanti-reads-integration.md) |
| Diagnose an error or empty/failed hub | [`references/troubleshooting.md`](references/troubleshooting.md) |

Runnable example invocations (reference + custom-genome) and the palette JSON format live
in [`assets/`](assets/).

## Conventions & footguns (high-value — internalize these)

- **UCSC command-line tools must be on PATH.** Always required: `gtfToGenePred`,
  `genePredToBed`, `bedToBigBed`. `twoBitInfo` is required with `--twobit`. `ixIxx` is
  needed for the Trix search index — if it's absent, Trix is **silently skipped** (no
  error, just no search). Install via `environment.yml` (bioconda) or `install_ucsc_tools.sh`.
- **Validate before uploading:** `hubCheck -noTracks example_output/hub.txt`.
- **Host the *raw* `hub.txt` URL** (GitHub: the `raw.githubusercontent.com` URL, not the
  repo page). Then UCSC → My Data → Track Hubs → Connected Hubs → paste the URL.
- **Inputs beyond the genome are dropped, with warnings.** CAGE/polyA/STAR/reference and
  chrom.sizes are filtered to chromosomes and coordinate bounds present in the genome — out
  of range features are silently removed.
- **Default UCSC filter caps:** exon count ≤ 400 and transcript length ≤ 150 kb. Larger
  genes (e.g. TTN) may hit these range limits.
- **HTML tables (`--tables`) need internet** to render — jQuery/DataTables load from a CDN.
- **Debug flags:** `--validate-only` (check tools+inputs and exit), `--dry-run` (build
  intermediates and stop before bigBed/hub), `--keep-temp` (keep temp dir; it's also kept
  automatically on failure).
- **`--install-only` is a test-harness flag, not a CLI flag.** It belongs to
  `python tests/test_sqanti_browser.py --install-only` (quick env check). The
  `sqanti_browser` command has no `--install-only` — for a main-CLI environment check use
  `--validate-only`.
- **Multiple hubs / samples:** give each run a distinct `--hub-name` so track labels are
  prefixed and comparable when several hubs are loaded together.

## Common options (quick reference; full detail in cli-reference.md)

`--tables` interactive HTML tables · `--sort-by` order isoforms (default `basic`) ·
`--category-tracks FSM,ISM,NIC` limit per-category tracks · `--no-category-tracks` main
track only · `--hub-name NAME` label/prefix · `--my-palette FILE` colors · `--no-highlight`
· `--trix`/`--no-trix` search index · `--star-sj` / `--CAGE-peak` / `--polyA-peak` /
`--refGTF` validation & comparison tracks · `--twobit` non-reference genome ·
`--chrom-sizes` explicit sizes · `--validate-only` / `--dry-run` / `--keep-temp` debug.
