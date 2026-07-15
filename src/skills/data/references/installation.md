# Installation

SQANTI-browser needs Python (>=3.8, dep: `pandas>=1.3.0`) plus a set of UCSC command-line
tools on PATH. Pick one of the routes below.

## Route A — Conda (recommended, ships UCSC tools)

```bash
conda env create -f environment.yml
conda activate sqanti-browser
pip install -e .          # install the package itself (env file does NOT install it)
```

`environment.yml` pulls `python>=3.8`, `pandas>=1.3.0`, and the UCSC bioconda packages
(`ucsc-bedtobigbed`, `ucsc-gtftogenepred`, `ucsc-genepredtobed`, `ucsc-hubcheck`,
`ucsc-twobitinfo`, `ucsc-fatotwobit`, `ucsc-ixixx`). If `ucsc-ixixx` fails to resolve,
install it via `install_ucsc_tools.sh` or skip Trix (search index is optional).

## Route B — pip + tool installer (no conda)

```bash
bash install_ucsc_tools.sh   # downloads UCSC binaries to /usr/local/bin (or ~/bin if unwritable)
pip install -e .             # or: pip install -r requirements.txt  (deps only)
```

`install_ucsc_tools.sh` auto-detects OS/arch and fetches: gtfToGenePred, genePredToBed,
bedToBigBed, hubCheck, bigBedInfo, ixIxx, twoBitInfo. It does NOT fetch `faToTwoBit` —
grab that separately if you need to build `.2bit` files (`conda install -c bioconda
ucsc-fatotwobit`, or the UCSC download page). If it installs to `~/bin`, add it to PATH:
`export PATH="$HOME/bin:$PATH"`.

## Route C — containers

Docker image on Docker Hub (v1.1.0+): `docker pull anaconesalab/sqanti-browser`. Or build:

```bash
docker build -t sqanti-browser:latest .
docker run -v $(pwd):/data sqanti-browser:latest \
    --gtf /data/corrected.gtf --classification /data/classification.txt \
    --output /data/my_hub --genome hg38
```

Singularity (`Singularity.def` in repo root):

```bash
singularity build sqanti-browser.sif Singularity.def
singularity run sqanti-browser.sif --gtf corrected.gtf \
    --classification classification.txt --output my_hub --genome hg38
```

Both containers build the conda env and set the entrypoint to `sqanti_browser`, so pass
flags directly.

## Which UCSC tools are required

| Tool | When needed |
|------|-------------|
| gtfToGenePred, genePredToBed, bedToBigBed | Always |
| twoBitInfo | Only with `--twobit` |
| ixIxx | Only for Trix search index — silently skipped if absent |
| faToTwoBit | Only to create a `.2bit` from FASTA (not at runtime) |
| hubCheck | Recommended, for validating the finished hub |

## Three ways to run

| Command | When |
|---------|------|
| `python -m sqanti_browser` | Recommended — uses the active (e.g. conda) Python, avoiding wrong-Python/arch issues. |
| `sqanti_browser` | After `pip install -e .` with the correct Python. |
| `python sqanti_browser.py` | From repo root; set `PYTHONPATH=.` if you hit `ModuleNotFoundError: No module named 'src'`. |

## Verify the environment

```bash
python -m sqanti_browser --validate-only \
    --gtf example/SQANTI3_QC_output/example_corrected.gtf \
    --classification example/SQANTI3_QC_output/example_classification.txt \
    --output /tmp/check --genome hg38
```

`--validate-only` checks tools + inputs and exits. (Note: `--install-only` is a flag of the
test harness `python tests/test_sqanti_browser.py --install-only`, not of the `sqanti_browser`
CLI — for a CLI env check use `--validate-only`.)

See [cli-reference.md](cli-reference.md) for the full flag list.
