# SQANTI-browser

SQANTI-browser is a tool that converts SQANTI3 transcriptome analysis outputs into UCSC Genome Browser track hubs for interactive visualization.

## Features

- 🎨 **Color-coded transcripts** by structural category
- 🔍 **Advanced filtering** with dropdowns and range sliders
- 📊 **Per-category tracks** for easy exploration
- 🔎 **Trix search** for finding isoforms using keywords for any attribute
- 📋 **Interactive HTML tables** with export capabilities
- 🧬 **SQANTI-reads compatible** - Process multi-sample experiments

## Documentation

📖 **[See the Wiki](../../wiki/home)** for detailed documentation:

**Quick Start:**
- 🚀 **[Quick Reference](../../wiki/quick_reference)** - One-page cheat sheet
- ❓ **[FAQ](../../wiki/faq)** - Common questions answered
- [Installation Guide](../../wiki/installation)
- [Usage Examples](../../wiki/usage_examples)
- [Hosting Guide](../../wiki/hosting) - Upload your hub

**Advanced Features:**
- [Interactive HTML Tables](../../wiki/html_tables) - Offline data exploration
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

## 🔄 Workflow

```
SQANTI3 Output → SQANTI-browser → web-access → UCSC Browser
```

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
2. Go to [UCSC Genome Browser](https://genome.ucsc.edu/) → **My Data** → **Track Hubs** → **Connected Hubs**
3. Enter the URL to your publicly available `hub.txt` file into the URL window
4. Click **Add Hub**

## Common Options

| Option | Description |
|--------|-------------|
| `--tables` | Generate interactive HTML tables for each category |
| `--sort-by none` | Sort isoforms by expression (default: `none`). Options: `length`, `FL`, `iso_exp` etc. |
| `--no-category-tracks` | Only generate the main track |
| `--star-sj SJ.out.tab` | Include STAR splice junctions track |
| `--CAGE_peak` | Include CAGE peaks for TSS validation (requires BED file) |
| `--polyA_peak` | Include polyA peaks for TTS validation (requires BED file) |
| `--refGTF` | Include reference annotation for direct comparison (requires GTF file) |

> 💡 **Tip:** For information on where to obtain and how to format CAGE and polyA peak files, see the [SQANTI3 wiki](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control#incorporating-cage-peak-data---cage_peak). See also our [Command Line Reference](wiki/command_line_reference.md#validation-tracks) for details.

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

