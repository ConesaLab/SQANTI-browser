# Non-reference genomes (.2bit)

For assemblies not hosted at UCSC (plants, non-model organisms, de novo assemblies, spike-in
sequences), supply the genome sequence as a `.2bit` file via `--twobit`. Requires the
`twoBitInfo` UCSC tool on PATH.

```bash
python -m sqanti_browser \
    --gtf corrected.gtf --classification classification.txt \
    --output my_hub --genome my_species_v1 \
    --twobit genome.2bit
```

`--genome` becomes the assembly name / genome subdir; it does not need to match any UCSC
assembly. No `--chrom-sizes` needed — chromosome sizes are extracted from the `.2bit`.

## What `--twobit` does automatically

1. Computes chrom.sizes from the `.2bit` via `twoBitInfo`.
2. Copies the `.2bit` into the hub at `<output>/<genome>/<genome>.2bit` (portable, self-
   contained).
3. Writes `genomes.txt` with `twoBitPath`, `organism`/`description`/`scientificName` set to
   user-defined, a `defaultPos` (first ~5 kb of the largest chromosome), and `htmlPath`.
4. Generates a genome description page `<genome>_genome.html` listing chromosomes and lengths.

## Creating a .2bit from FASTA

```bash
faToTwoBit genome.fasta genome.2bit
```

`faToTwoBit` is not fetched by `install_ucsc_tools.sh` — install it with
`conda install -c bioconda ucsc-fatotwobit` (it is in the conda `environment.yml`).

## Worked example — SIRVs

The repo ships a spike-in example under `example/SQANTI3_QC_custom_genome/` with
`SIRVS.2bit`, `sirv_corrected.gtf`, `sirv_classification.txt`:

```bash
python -m sqanti_browser \
    --gtf example/SQANTI3_QC_custom_genome/sirv_corrected.gtf \
    --classification example/SQANTI3_QC_custom_genome/sirv_classification.txt \
    --output sirv_hub --genome SIRV \
    --twobit example/SQANTI3_QC_custom_genome/SIRVS.2bit \
    --tables
```

Then host and load like any hub (see [hosting.md](hosting.md)).

## Troubleshooting

| Problem | Fix |
|---------|-----|
| `twoBitInfo not found` | Install UCSC tools (`bash install_ucsc_tools.sh` or conda). |
| Invalid `.2bit` | Recreate from FASTA with `faToTwoBit`. |
| Coords dropped | Validation/reference features off the genome's chromosomes are filtered out (expected). |
