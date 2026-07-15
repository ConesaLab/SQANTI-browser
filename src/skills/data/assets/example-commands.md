# Ready-to-run example invocations

These use the data bundled in the SQANTI-browser repo (`example/`). Run from the repo root.
Both require the UCSC command-line tools on PATH (see `references/installation.md`).

## Reference genome (hg38) — the canonical smoke test

```bash
python -m sqanti_browser \
  --gtf example/SQANTI3_QC_output/example_corrected.gtf \
  --classification example/SQANTI3_QC_output/example_classification.txt \
  --output example_output \
  --genome hg38
```

With the full set of validation tracks + interactive tables:

```bash
python -m sqanti_browser \
  --gtf example/SQANTI3_QC_output/example_corrected.gtf \
  --classification example/SQANTI3_QC_output/example_classification.txt \
  --output example_output \
  --genome hg38 \
  --refGTF example/SQANTI3_QC_output/gencode.v38.basic_chr22.gtf \
  --star-sj example/SQANTI3_QC_output/exampleSJ.out.tab \
  --CAGE-peak example/SQANTI3_QC_output/chr22.human.refTSS_v3.1.hg38.bed \
  --polyA-peak example/SQANTI3_QC_output/polyApeaks.atlas.GRCh38.bed \
  --tables --sort-by iso_exp
```

## Non-reference genome (SIRVs, via .2bit)

```bash
python -m sqanti_browser \
  --gtf example/SQANTI3_QC_custom_genome/sirv_corrected.gtf \
  --classification example/SQANTI3_QC_custom_genome/sirv_classification.txt \
  --output example_output/custom_genome \
  --genome SIRVs \
  --twobit example/SQANTI3_QC_custom_genome/SIRVS.2bit
```

## Validate the hub before uploading

```bash
hubCheck -noTracks example_output/hub.txt
```

## Custom palette

Copy `example/example_palette.json` (top-level keys `standard`, `highlight`,
`validation_tracks`; RGB arrays or hex) and pass it with `--my-palette`. Partial palettes
are allowed — unspecified categories fall back to defaults. See
`references/custom-coloring.md`.
