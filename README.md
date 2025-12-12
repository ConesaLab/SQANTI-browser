# SQANTI-browser

Convert SQANTI3 transcriptome analysis outputs into UCSC Genome Browser track hubs for interactive visualization.

## Features

- 🎨 **Color-coded transcripts** by structural category
- 🔍 **Advanced filtering** with dropdowns and range sliders
- 📊 **Per-category tracks** for easy exploration
- 🔎 **Trix search** for finding isoforms by any attribute
- 📋 **Interactive HTML tables** with export capabilities

## Quick Start

### 1. Install Dependencies

```bash
# Install UCSC tools
bash install_ucsc_tools.sh

# Install Python dependencies
pip install -r requirements.txt
```

### 2. Run the Converter

```bash
python sqanti3_to_UCSC.py \
    --gtf your_corrected.gtf \
    --classification your_classification.txt \
    --output my_hub \
    --genome hg38
```

### 3. Upload to UCSC

1. Upload the output directory to a web server (e.g., GitHub)
2. Go to [UCSC Genome Browser](https://genome.ucsc.edu/) → **My Data** → **Track Hubs**
3. Enter the URL to your `hub.txt` file
4. Click **Add Hub**

## Common Options

| Option | Description |
|--------|-------------|
| `--tables` | Generate interactive HTML tables for each category |
| `--github-repo user/repo` | Auto-generate GitHub raw URLs |
| `--sort-by iso_exp` | Sort isoforms by expression (or `length`, `FL`, etc.) |
| `--no-category-tracks` | Only generate the main track |
| `--star-sj SJ.out.tab` | Include STAR splice junctions track |

## Color Legend

| Category | Color | Hex Code |
|----------|-------|----------|
| Full-splice match (FSM) | Blue | <span style="color:#6BAED6">#6BAED6</span> |
| Incomplete-splice match (ISM) | Orange | <span style="color:#FC8D59">#FC8D59</span> |
| Novel in catalog (NIC) | Green | <span style="color:#78C679">#78C679</span> |
| Novel not in catalog (NNC) | Red | <span style="color:#EE6A50">#EE6A50</span> |
| Genic | Gray | <span style="color:#969696">#969696</span> |
| Antisense | Teal | <span style="color:#66C2A4">#66C2A4</span> |
| Fusion | Gold | <span style="color:#DAA520">#DAA520</span> |
| Intergenic | Salmon | <span style="color:#E9967A">#E9967A</span> |
| Genic intron | Cyan | <span style="color:#41B6C4">#41B6C4</span> |

## Documentation

📖 **[See the Wiki](../../wiki)** for detailed documentation:

- [Installation Guide](../../wiki/Installation)
- [Usage Examples](../../wiki/Usage-Examples)
- [Filtering in UCSC](../../wiki/Filtering-in-UCSC)
- [Trix Search Syntax](../../wiki/Trix-Search)
- [Isoform Ordering](../../wiki/Isoform-Ordering)
- [Command Line Reference](../../wiki/Command-Line-Reference)
- [Output File Formats](../../wiki/Output-Files)
- [Troubleshooting](../../wiki/Troubleshooting)

## License

GNU General Public License v3.0 - see [LICENSE](LICENSE)

## Citation

If you use this tool, please cite:
- [SQANTI3](https://github.com/ConesaLab/SQANTI3)
- [UCSC Genome Browser](https://genome.ucsc.edu/)
