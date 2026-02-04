<img width="952" height="309" alt="image" src="https://github.com/user-attachments/assets/df680d97-ed84-4657-972e-25fb61471c57" />
  
  
**SQANTI-browser** is a tool that converts SQANTI3 transcriptome analysis outputs into UCSC Genome Browser track hubs for interactive visualization.

## Features

- 🎨 **Color-coded transcripts** by structural category (custom palettes via `--my-palette`)
- 🔍 **Advanced filtering** with dropdowns and range definitions in UCSC
- 📊 **Per-category tracks** for easy exploration
- 🔎 **Trix search** for finding isoforms by attribute
- 📋 **Interactive HTML tables** with multi-select rows, export, **Trix string generation**, and **Generate Filter String** (newline-separated list for Table Browser)
- 👁 **View specific isoforms** via Table Browser (Custom Tracks)
- 🎯 **Validation tracks** (CAGE, PolyA, STAR, reference)
- 🧬 **SQANTI-reads compatible** – Process multi-sample experiments

## Documentation

📖 **[See the Wiki](../../wiki/home)** for detailed documentation.

All docs live in `wiki/`. Edit there, review, then push when ready (assistants never push to GitHub). See `wiki/README.md` for the workflow.

**Quick Start:**
- 🚀 **[Quick Reference](../../wiki/quick_reference)** - One-page cheat sheet
- ❓ **[FAQ](../../wiki/faq)** - Common questions answered
- [Installation Guide](../../wiki/installation)
- [Usage Examples](../../wiki/usage_examples)
- [Hosting Guide](../../wiki/hosting) - Upload your hub

**Advanced Features:**
- [Interactive HTML Tables](../../wiki/html_tables) - Offline data exploration
- [Custom Color Palettes](../../wiki/custom_palette) - Palette format, defaults, and validation track colors
- [Working with Custom Genomes](../../wiki/custom_genomes) (2bit files)
- [SQANTI-reads Integration](../../wiki/sqanti_reads_integration) - Multi-sample experiments
- [Filtering in UCSC](../../wiki/filtering_in_UCSC)
- [Trix Search Syntax](../../wiki/trix_search)
- [Isoform Ordering](../../wiki/isoform_ordering)

**Reference:**
- [Command Line Reference](../../wiki/command_line_reference)
- [Output File Formats](../../wiki/output_files)
- [Glossary](../../wiki/glossary)
- [Troubleshooting](../../wiki/troubleshooting)

## Project Structure

```
SQANTI-browser/
├── sqanti_browser.py       # Main converter script
├── install_ucsc_tools.sh   # UCSC tools installer
├── src/
│   ├── bed_processor.py    # BED conversion logic
│   ├── filter_isoforms.py  # Interactive HTML table/report generator
│   ├── hub_generator.py    # UCSC hub file generation
│   ├── validation_tracks.py# CAGE, PolyA, STAR, reference tracks
│   ├── constants.py        # Filter limits, palette defaults
│   └── utils.py            # Shared utilities
├── example/
│   ├── example_usage.py    # Example workflow script
│   ├── example_palette.json
│   └── SQANTI3_QC_output/  # Example data
├── tests/
│   ├── test_sqanti_browser.py   # Integration tests
│   ├── test_unit.py             # Unit tests
│   └── test_edge_cases.py       # Edge case tests
└── example_output/         # Sample hub output (pre-generated)
```

## 🔄 Workflow

```
SQANTI3 Output → SQANTI-browser → web-access → UCSC Browser
```

## Quick Start

### 1. Install Dependencies

```bash
# Install UCSC tools
bash install_ucsc_tools.sh

# Install Python package (recommended)
pip install -e .

# Or install dependencies only
pip install -r requirements.txt
```

**How to run:** You can use any of these (run from project root unless you used `pip install -e .`):

- `python -m sqanti_browser` — Recommended. Uses whatever `python` is active (e.g. conda), so you avoid wrong-Python / architecture issues.
- `sqanti_browser` — After `pip install -e .`, if you installed with the correct Python.
- `python sqanti_browser.py` — Uses the `python` in your PATH; needs project root in `PYTHONPATH` or package installed.

### 2. Run the Converter

```bash
python -m sqanti_browser \
    --gtf your_corrected.gtf \
    --classification your_classification.txt \
    --output my_hub \
    --genome hg38
```

Run the example workflow with bundled data:

```bash
python example/example_usage.py
```

### 3. Upload to UCSC

1. Upload the output directory to a web server (e.g., GitHub)
2. Go to [UCSC Genome Browser](https://genome.ucsc.edu/) → **My Data** → **Track Hubs** → **Connected Hubs**
3. Enter the URL to your publicly available `hub.txt` file into the URL window
4. Click **Add Hub**

## Common Options

| Option | Description |
|--------|-------------|
| `--tables` | Generate interactive HTML tables for each category |
| `--sort-by none` | Sort isoforms (Default `none`). Options: `FL`, `iso_exp`, `length`, etc. |
| `--no-category-tracks` | Only generate the main track |
| `--category-tracks FSM,ISM,NIC` | Only create tracks for these categories (abbreviated: FSM, ISM, NIC, NNC, antisense, genic_intron, genic_genomic, intergenic, fusion) |
| `--no-highlight` | Disable highlight coloring for top FL isoforms |
| `--my-palette FILE` | Custom color palette JSON file (see [Custom palette wiki](../../wiki/custom_palette)) |
| `--star-sj SJ.out.tab` | Include STAR splice junctions track |
| `--CAGE-peak` | Include CAGE peaks for TSS validation (requires BED file) |
| `--polyA-peak` | Include polyA peaks for TTS validation (requires BED file) |
| `--refGTF` | Include reference annotation for direct comparison (requires GTF file) |

### Testing

```bash
# Full integration tests (requires UCSC tools)
python tests/test_sqanti_browser.py

# Unit and edge case tests (no UCSC tools needed)
python tests/test_unit.py -v
python tests/test_edge_cases.py -v
```

Use `--install-only` for a quick environment check (UCSC tools, Python deps).

### Filter Reports (standalone)

Generate interactive HTML reports from a classification file:

```bash
python src/filter_isoforms.py --classification your_classification.txt --output-dir report_output
```

> 💡 **Tip:** For information on where to obtain and how to format CAGE and polyA peak files, see the [SQANTI3 wiki](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control#incorporating-cage-peak-data---cage_peak). See also our [Command Line Reference](../../wiki/command_line_reference) for details.

## Viewing Specific Isoforms

To view only a subset of isoforms (e.g., a few of interest): use the **Interactive HTML tables** (`--tables`), select the rows you want (click to select multiple), click **Generate Filter String**, then in UCSC go to **Tools → Table Browser**, choose your hub’s SQANTI3 track, use **Identifiers → paste list** for the newline-separated IDs, set the output format to **custom track**, and load it in the Genome Browser. See the [Filtering in UCSC](../../wiki/filtering_in_UCSC) wiki page for details.

## Trix Search

The hub includes a Trix search index for finding isoforms by attribute. **All search terms use underscores** (e.g., `structural_category_full_splice_match`, `strand_plus`, `coding_coding`). Hyphens and spaces in category names are automatically normalized, so you can search with `full_splice_match` instead of remembering `full-splice_match` (although this will also work). Use the "Generate Trix String" button in the interactive HTML tables to get easy ready-to-paste search terms.

## Default Track Display

By default, the **main SQANTI3 track** opens visible in **pack mode**. **Validation tracks** (STAR junctions, CAGE peaks, PolyA peaks, reference)—when included—also show by default. **Category-specific tracks** (e.g., full-splice_match, novel_in_catalog) start **hidden** so the view stays uncluttered. You can turn on individual category tracks from the track list as needed. Track description pages use left-aligned, transparent styling to match the UCSC Genome Browser interface.

For default color names (e.g. steelblue, coral, peru), custom palette format, validation track colors, and examples, see the **[Custom Color Palettes](../../wiki/custom_palette)** wiki page.

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for setup and guidelines.

## License

GNU General Public License v3.0 - see [LICENSE](LICENSE)

## Performance

| Dataset Size | Time | Memory |
|--------------|------|--------|
| ~1K transcripts | 1-2 min | 500MB |
| ~10K transcripts | 5-10 min | 1GB |
| ~100K transcripts | 30-60 min | 4GB |

## Citation

If you use this tool, please cite:
- https://github.com/conesalab/SQANTI-browser
- Perez, G. et al. **The UCSC Genome Browser database: 2025 update**. *Nucleic Acids Res* (2025). https://doi.org/10.1093/nar/gkae974
- Pardo-Palacios, F.J., Arzalluz-Luque, A. et al. **SQANTI3: curation of long-read transcriptomes for accurate identification of known and novel isoforms**. *Nat Methods* (2024). https://doi.org/10.1038/s41592-024-02229-2

If you used SQANTI-reads for multi-sample visualization, please cite:
- Keil N, Monzó C, McIntyre L, Conesa A (2025). **Quality assessment of long read data in multisample lrRNA-seq experiments with SQANTI-reads**. *Genome Res*, 35(4), 987. https://doi.org/10.1101/gr.280021.124

