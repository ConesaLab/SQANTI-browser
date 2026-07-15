# Inputs & outputs

## Inputs

**Required** (both from a SQANTI3 QC run):

- `--gtf` — SQANTI3 **corrected GTF** (`*_corrected.gtf`): transcript models.
- `--classification` — SQANTI3 **classification table** (`*_classification.txt`, tab-
  separated). The `isoform` column is used as the feature `name`; ~48 columns (e.g.
  `structural_category`, `associated_gene`, `associated_transcript`, `FL`, `iso_exp`,
  `coding`, `diff_to_TSS`, `dist_to_CAGE_peak`, …) are carried into the bigBed as extra
  fields and become UCSC filters + Trix search tokens.

**Optional:**

- `--chrom-sizes` (tab file) or `--twobit` (`.2bit`) — genome sizes / non-reference genome.
- `--star-sj` (STAR `SJ.out.tab`), `--CAGE-peak` (BED), `--polyA-peak` (BED),
  `--refGTF` (reference GTF) — validation/comparison tracks.

Everything except the transcripts is clipped to the genome's chromosomes and bounds; out-
of-range rows are dropped with warnings. (v1.1.1+ derives chrom.sizes considering
validation inputs so `bedToBigBed` won't fail when junction/peak coords exceed transcript
extents.)

## Output — a UCSC track hub

```
<output>/
├── hub.txt                 # hub entry point (name, email, references)  <- load THIS url in UCSC
├── genomes.txt             # assembly → trackDb mapping (+ custom genome def if --twobit)
├── README.md               # hub docs with a color legend
└── <genome>/
    ├── <genome>_sqanti3.bb           # main track: all transcripts, colored by category (bigBed 12+N)
    ├── <genome>_sqanti3_<category>.bb# one per structural category (hidden by default)
    ├── <genome>_star_sj.bb           # optional, with --star-sj
    ├── <genome>_cage_peaks.bb        # optional, with --CAGE-peak
    ├── <genome>_polya_peaks.bb       # optional, with --polyA-peak
    ├── <genome>_reference.bb         # optional, with --refGTF
    ├── trackDb.txt                   # track defs, filters, display modes
    ├── groups.txt                    # groups: transcripts (1), validation (2)
    ├── trix.ix / trix.ixx            # search index (omitted if ixIxx missing / --no-trix)
    ├── <genome>_sqanti3_track.html   # per-track description pages
    ├── <genome>.2bit                 # only with --twobit
    └── <genome>_genome.html          # only with --twobit
```

With `--tables`, also `<output>/table_reports/`:
`complete_transcriptome_isoforms.html` + one `<category>_isoforms.html` per category
(DataTables-based; needs internet for CDN jQuery/DataTables).

## Track display defaults

- Main `<genome>_sqanti3` track: **visible, pack mode**, `itemRgb on`, searchable
  (`searchIndex name`, `searchTrix trix.ix`).
- Per-category tracks: **hidden** (turn on as needed).
- Validation tracks (STAR/CAGE/polyA/reference): shown (full/dense) when included.

## bigBed format

BED12 + N extra columns. Base BED12: `chrom, chromStart, chromEnd, name, score, strand,
thickStart, thickEnd, itemRgb, blockCount, blockSizes, chromStarts` (thick = CDS,
itemRgb = category color). N extra columns = the classification fields, exposed as UCSC
filters and Trix tokens.

## Rough sizing

Hub ≈ 50% of input size; `.bb` ≈ 10–50% of the GTF; trix ≈ 1–5% of the classification.
Performance: ~1K transcripts 1–2 min/500 MB; ~100K transcripts 30–60 min/4 GB.
