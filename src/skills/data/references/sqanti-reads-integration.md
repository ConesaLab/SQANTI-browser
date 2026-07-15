# SQANTI-reads integration

SQANTI-reads (part of the SQANTI3 toolkit) does multi-sample QC of long-read RNA-seq. Its
per-sample outputs feed straight into SQANTI-browser. **Build one hub per sample, not a
merged hub** — each sample keeps independent filtering and can be shared separately, and you
load several hubs together in UCSC to compare.

## Inputs per sample

- **GTF**: the SQANTI3 corrected GTF (`<sample>_corrected.gtf`), unchanged by SQANTI-reads.
- **Classification**: `<sample>_reads_classification.txt` — same structure as a normal
  SQANTI3 classification file plus two extra columns, `jxn_string` (junction chain) and
  `jxnHash`. SQANTI-browser handles it identically and carries both extra columns into the
  bigBed, where they become searchable/filterable.

## Use a distinct `--hub-name` per sample

Without `--hub-name`, every hub uses the same track labels (e.g. "SQANTI3 Transcripts"), so
loading multiple samples shows identical names. Give each sample a distinct `--hub-name` so
every track label is prefixed and comparable in the track list:

```bash
python -m sqanti_browser \
    --gtf SAMPLE001/sample1_PacBio_corrected.gtf \
    --classification SAMPLE001/sample1_PacBio_reads_classification.txt \
    --output hubs/sample1_PacBio_hub \
    --genome hg38 --hub-name "sample1_PacBio" --tables
```

`--hub-name "sample1_PacBio"` yields labels like "sample1_PacBio SQANTI3 Transcripts".

## Batch processing

Loop over a design CSV, reusing a common reference GTF via `--refGTF` (the same annotation
works for all same-genome samples):

```bash
DESIGN_CSV="design.csv"; INPUT_DIR="sqanti_reads_output"
OUTPUT_BASE="hubs"; GENOME="hg38"; REF_GTF="gencode.v38.annotation.gtf"

tail -n +2 "$DESIGN_CSV" | while IFS=',' read -r sampleID file_acc rest; do
    GTF="${INPUT_DIR}/${file_acc}/${sampleID}_corrected.gtf"
    CLS="${INPUT_DIR}/${file_acc}/${sampleID}_reads_classification.txt"
    if [[ -f "$GTF" && -f "$CLS" ]]; then
        python -m sqanti_browser \
            --gtf "$GTF" --classification "$CLS" \
            --output "${OUTPUT_BASE}/${sampleID}_hub" \
            --genome "$GENOME" --hub-name "$sampleID" \
            --refGTF "$REF_GTF" --tables
    fi
done
```

Add `--star-sj <file>` per sample if STAR `SJ.out.tab` files are available.

## Hosting and comparison

Each hub is a self-contained directory (`hub.txt`, `genomes.txt`, `<genome>/`). Upload each
to a public location and add its raw `hub.txt` URL in UCSC (see [hosting.md](hosting.md)):

```
https://raw.githubusercontent.com/USER/REPO/main/sample1_PacBio_hub/hub.txt
https://raw.githubusercontent.com/USER/REPO/main/sample2_PacBio_hub/hub.txt
```

All samples then appear in the track list, prefixed by `--hub-name`, for side-by-side
comparison, independent filtering, and figure export.

## Notes

- Different genomes across samples → process separately with the right `--genome`
  (`hg38`, `mm10`, …); UCSC handles multiple assemblies.
- Merging into one hub is possible only by hand-editing `trackDb.txt`; separate hubs are
  recommended.
