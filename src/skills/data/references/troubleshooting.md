# Troubleshooting

## Quick checks

```bash
python -m sqanti_browser ... --validate-only   # check tools + inputs, then exit
hubCheck -noTracks hub.txt                      # validate generated hub structure
```

(The `sqanti_browser` CLI has no `--install-only`; that flag is on the test harness
`python tests/test_sqanti_browser.py --install-only`. For a CLI check use `--validate-only`.)

## Missing UCSC tools

`gtfToGenePred: command not found` (or genePredToBed / bedToBigBed / twoBitInfo). Install
and put on PATH:

```bash
bash install_ucsc_tools.sh        # or: conda env create -f environment.yml
export PATH="$HOME/bin:$PATH"      # if installed to ~/bin
```

gtfToGenePred, genePredToBed, bedToBigBed are always required; twoBitInfo only with
`--twobit`.

## Wrong Python / missing pandas

`ModuleNotFoundError: No module named 'pandas'` (or `'src'`) usually means the wrong
interpreter. Prefer `python -m sqanti_browser` so the active (conda) Python is used. Install
deps with `pip install -r requirements.txt`. For `No module named 'src'` when running
`python sqanti_browser.py`, run from the repo root with `PYTHONPATH=.`.

## Trix search skipped

`ixIxx` is optional and **silently skipped** if not on PATH — the hub builds but has no
search index. If search returns nothing: confirm `trix.ix` and `trix.ixx` exist in the
`<genome>/` dir, that `trackDb.txt` has `searchTrix`, install `ixIxx`, then regenerate.

## Coordinates outside the genome dropped

`Filtered CAGE peaks: 100/10000 peaks kept` and similar are **normal**. Validation/reference
features on chromosomes not present in your transcriptome (or beyond chromosome bounds) are
dropped. If your data is chr22 only, other chromosomes' features are filtered out.

## bedToBigBed fails / "file not sorted"

Check chromosome names match the assembly, remove special characters from transcript IDs, and
inspect intermediates with `--keep-temp`.

## Filter ranges look capped

Track range sliders default to the caps in `src/constants.py` (e.g. `exons` 0–400, `length`
0–150,000). Values beyond the default range aren't reachable by the slider — widen the box
manually. See [filtering-and-search.md](filtering-and-search.md).

## Filters not appearing / not applying

Right-click the track → **Configure** (not plain settings); set visibility to full/pack (not
hide); click **Submit** after setting filters; clear cache by adding `udcTimeout=5` to the
URL and reconnecting.

## Tables don't render

The HTML tables load DataTables/jQuery from a **CDN** — the browser needs internet access
when opening them.

## hubCheck errors

| Error | Fix |
|-------|-----|
| `can't open hub.txt` | Check path/permissions; for remote, confirm the URL is public. |
| `genome not found` | `genomes.txt` must reference the correct assembly. |
| `trackDb.txt missing` | Must exist in the `<genome>/` dir. |
| `bigDataUrl not found` / `can't open bigBed` | The `.bb` file is missing at its path — re-upload or regenerate. |

## Debug flags

| Flag | Effect |
|------|--------|
| `--validate-only` | Validate tools + inputs, exit before building. |
| `--dry-run` | Build the intermediate merged BED, exit before bigBed/hub. |
| `--keep-temp` | Preserve the temp dir (also kept automatically on failure) for inspection. |

See [cli-reference.md](cli-reference.md) and [installation.md](installation.md).
