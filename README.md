# SQANTI-browser

SQANTI-browser converts SQANTI3 transcriptome analysis outputs into UCSC Genome Browser track hubs. It takes SQANTI3 `*_corrected.gtf` and `*_classification.txt`, produces bigBed files, and generates all hub configuration needed for immediate visualization in the UCSC Genome Browser.

## Features

- **Automatic Conversion**: Converts SQANTI3 GTF and classification files to bigBed format
- **Color Coding**: Transcripts are automatically colored by structural category using the same color scheme as SQANTI3
- **Per-Category Tracks**: Separate bigBed tracks for each structural category for easy filtering
- **Trix Search Index**: Full-text search for isoforms by any classification attribute (gene, category, strand, etc.)
- **Interactive HTML Reports**: Filterable tables for each structural category with range filters, dropdowns, and Trix string generation
- **Hub Generation**: Creates complete UCSC Genome Browser hub files with tracks visible by default
- **Comprehensive Logging**: Detailed progress and error reporting
- **Advanced**: Is compatible with user-defined chrom.sizes, .2bit, can also format and include STAR junction files.

## Prerequisites

### Required Software

1. **UCSC Tools** (required for file format conversion and search indexing):
   - `gtfToGenePred`
   - `genePredToBed` 
   - `bedToBigBed`
   - `ixIxx`

   **Installation Options:**

   **Option A: Use the installation script for automatic installation on Linux/Mac**
   ```bash
   sh install_ucsc_tools.sh
   ```
   
   **Option B: Download from UCSC**
   ```bash
   # Download from UCSC (Linux/Mac)
   wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
   wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
   wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
   wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ixIxx
   
   # Make executable
   chmod +x gtfToGenePred genePredToBed bedToBigBed ixIxx
   
   # Add to PATH or move to a directory in your PATH
   sudo mv gtfToGenePred genePredToBed bedToBigBed ixIxx /usr/local/bin/
   ```
   
   **Option C: Conda installation**
   ```bash
   conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed ucsc-bedtobigbed ucsc-ixixx
   ```

2. **Python Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Basic Usage

```bash
python sqanti3_to_UCSC.py \
    --gtf your_transcriptome_corrected.gtf \
    --classification your_transcriptome_classification.txt \
    --output output_directory \
    --genome hg38
```

### With Interactive HTML Tables

Generate both the track hub and interactive HTML reports for each structural category:

```bash
python sqanti3_to_UCSC.py \
    --gtf your_transcriptome_corrected.gtf \
    --classification your_transcriptome_classification.txt \
    --output output_directory \
    --genome hg38 \
    --tables
```

This creates a `table_reports/` subdirectory with filterable HTML tables for each structural category.

### Advanced Usage with Custom Chromosome Sizes

If you have a pre-existing chromosome sizes file:

```bash
python sqanti3_to_UCSC.py \
    --gtf your_transcriptome_corrected.gtf \
    --classification your_transcriptome_classification.txt \
    --output output_directory \
    --genome hg38 \
    --chrom-sizes /path/to/chrom.sizes
```
Alternatively, if you have a genome `.2bit` file, you can derive `chrom.sizes` using `twoBitInfo` automatically:

```bash
python sqanti3_to_UCSC.py \
    --gtf your_transcriptome_corrected.gtf \
    --classification your_transcriptome_classification.txt \
    --output output_directory \
    --genome hg38 \
    --twobit /path/to/genome.2bit
```
Priority order for chromosome sizes resolution:
- `--chrom-sizes` (highest)
- `--twobit` (via twoBitInfo)
- Extract from GTF (fallback)

### GitHub Integration (Recommended for UCSC Genome Browser)

For automatic GitHub raw URL generation, specify your repository:

```bash
python sqanti3_to_UCSC.py \
    --gtf your_transcriptome_corrected.gtf \
    --classification your_transcriptome_classification.txt \
    --output output_directory \
    --genome hg38 \
    --github-repo yourusername/your-repository-name
```

This generates hub files with GitHub raw URLs (`https://raw.githubusercontent.com/...`), ready to use directly in UCSC Genome Browser without manual editing.

### Optional Enhancements

- **STAR Splice Junctions Track**
  - Convert STAR `SJ.out.tab` into a bigBed junctions track and include it in the hub:
  ```bash
  python sqanti3_to_UCSC.py \
      --gtf your_corrected.gtf \
      --classification your_classification.txt \
      --output output_directory \
      --genome hg38 \
      --star-sj path/to/SJ.out.tab
  ```
  - This creates `{genome}_star_sj.bb` and adds a `STAR Junctions` track.

### Standalone HTML Reports

You can also generate HTML reports independently using `filter_isoforms.py`:

```bash
python filter_isoforms.py \
    --classification your_classification.txt \
    --output-dir output_reports
```

This generates one `.html` file per structural category (e.g., `full-splice_match_isoforms.html`, `novel_in_catalog_isoforms.html`).

**HTML Report Features:**
- **Interactive tables** with sorting and column-specific filtering
- **Range filtering** for numeric columns (e.g., type `100:1000` in the length column)
- **Dropdown filters** for categorical columns (strand, structural_category, coding, etc.)
- **Export to CSV/Excel**
- **Generate Trix String button** to create UCSC Genome Browser search queries
- **SVG schematics** showing each structural category visually
- **Documentation links** to SQANTI3 Wiki for column explanations
- Self-contained and shareable files

Use `--include-sequences` to keep the `ORF_seq` column (excluded by default).

### Validation and Dry Run

- **Validate only**: check tools and inputs, parse classification, then exit.
```bash
python sqanti3_to_UCSC.py \
    --gtf your_corrected.gtf \
    --classification your_classification.txt \
    --output output_directory \
    --genome hg38 \
    --validate-only
```

- **Dry run**: process up to the enhanced BED (with colors and encoded names), skip bigBed/hub creation.
```bash
python sqanti3_to_UCSC.py \
    --gtf your_corrected.gtf \
    --classification your_classification.txt \
    --output output_directory \
    --genome hg38 \
    --dry-run
```

- **Keep temporary files**: preserve intermediate files for debugging.
```bash
python sqanti3_to_UCSC.py \
    --gtf your_corrected.gtf \
    --classification your_classification.txt \
    --output output_directory \
    --genome hg38 \
    --keep-temp
```

### Command Line Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--gtf` | Yes | Path to SQANTI3 corrected GTF file |
| `--classification` | Yes | Path to SQANTI3 classification file |
| `--output` | Yes | Output directory for generated files (supports relative and absolute paths) |
| `--genome` | Yes | Genome assembly name (e.g., hg38, mm10, dm6) |
| `--chrom-sizes` | No | Optional path to chromosome sizes file |
| `--github-repo` | No | GitHub repository (format: username/repository) for automatic raw URL generation |
| `--github-branch` | No | GitHub branch (default: main) |
| `--star-sj` | No | Path to STAR `SJ.out.tab` to add a splice junctions track |
| `--tables` | No | Generate interactive HTML tables for each structural category |
| `--keep-temp` | No | Keep temporary files for debugging |
| `--validate-only` | No | Validate inputs and exit without processing |
| `--dry-run` | No | Process up to enhanced BED, skip bigBed/hub creation |

**Output Directory Flexibility:**
- **Relative paths**: `./my_project`, `../analysis`
- **Absolute paths**: `~/Documents/research`, `/home/username/projects`
- **Home directory**: `~/my_hub`, `$HOME/transcriptome_analysis`

**GitHub Integration:**
- **Automatic raw URLs**: When `--github-repo` is specified, all file references use GitHub raw URLs via `https://raw.githubusercontent.com/{repo}/{branch}/...`
- **No manual editing**: Hub files are ready to use directly in UCSC Genome Browser
- **Branch support**: Specify different branches with `--github-branch`

## Searching in UCSC Genome Browser

The tool generates a Trix search index that enables powerful text-based search in the UCSC Genome Browser.

### Search Syntax

Use the search box in the UCSC Genome Browser with the following formats:

| Search Term | Description | Example |
|-------------|-------------|---------|
| `isoform_id` | Search by exact isoform ID | `PB.1234.1` |
| `gene_name` | Search by associated gene | `BRCA1` or `ENSG00000012048` |
| `structural_category_X` | Filter by structural category | `structural_category_novel_in_catalog` |
| `strand_plus` / `strand_minus` | Filter by strand | `strand_plus` |
| `coding_coding` / `coding_non_coding` | Filter by coding status | `coding_coding` |
| `subcategory_X` | Filter by subcategory | `subcategory_mono-exon` |
| `bite_True` / `bite_False` | Filter by BITE status | `bite_True` |
| `predicted_NMD_True` | Filter by NMD prediction | `predicted_NMD_True` |
| `polyA_motif_found_True` | Filter by polyA motif | `polyA_motif_found_True` |

### Combining Search Terms

You can combine multiple terms with spaces (AND logic):

```
structural_category_novel_in_catalog strand_plus coding_coding
```

This finds all novel_in_catalog isoforms on the plus strand that are coding.

### Using the Generate Trix String Button

The HTML reports include a "Generate Trix String" button that helps create search queries:

1. **Filter-based**: Set your desired filters in the dropdowns, then click "Generate Trix String"
2. **Row-based**: Click on a specific isoform row to select it (highlighted in blue), then click "Generate Trix String" to get search terms for that isoform

The generated string can be copied directly into the UCSC Genome Browser search box.

> **Note:** Range filters (e.g., `100:1000` for length) are not supported by Trix search and will be ignored when generating the search string.

## Transcript Coloring

Transcripts are automatically colored based on their structural category:

| Structural Category | Color | Hex Code |
|---------------------|-------|----------|
| Full splice match | Blue | #6BAED6 |
| Incomplete splice match | Orange | #FC8D59 |
| Novel in catalog | Green | #78C679 |
| Novel not in catalog | Red | #EE6A50 |
| Genic | Gray | #969696 |
| Antisense | Teal | #66C2A4 |
| Fusion | Gold | #DAA520 |
| Intergenic | Salmon | #E9967A |
| Genic intron | Cyan | #41B6C4 |

## Uploading to UCSC Genome Browser

### Step 1: Host Your Files

Upload all generated files to a web-accessible location (e.g., GitHub Pages, your institution's web server, or cloud storage with public URLs).

### Step 2: Add the Hub

1. Go to the [UCSC Genome Browser](https://genome.ucsc.edu/)
2. Navigate to **My Data** → **Track Hubs**
3. Enter the URL to your `hub.txt` file
4. Click **Add Hub**

### Step 3: View Your Tracks

1. Choose the appropriate genome assembly from the dropdown
2. Your SQANTI3 tracks will appear in the track list (visible by default as "full")
3. Use the search box to find specific isoforms using Trix search syntax

## Example Workflow

### 1. Run SQANTI3

```bash
# Your existing SQANTI3 pipeline
python SQANTI3.py -g reference.gtf -j reference.gff3 -o output_dir input.fasta
```

### 2. Run the converter with HTML tables

```bash
python sqanti3_to_UCSC.py \
    --gtf output_dir/input_corrected.gtf \
    --classification output_dir/input_classification.txt \
    --output hub_output \
    --genome hg38 \
    --tables
```

### 3. Upload and View

1. Upload the `hub_output` directory to your web server (for example, your GitHub repository)
2. Check that the paths in your `hub.txt` file correspond to where your files are accessible in the web
3. Go to https://genome.ucsc.edu/index.html
4. Click on My Data -> Track Hubs -> My Hubs and add the link to your publicly available hub.txt file into the URL window.
4. Click on Add Hub.
5. Navigate to your region of interest or use the search box with Trix search syntax
6. View and filter transcripts by structural category

### Debug Mode

For detailed debugging, you can modify the logging level in the script:

```python
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
```

## File Format Details

### Input Files

- **GTF File**: Standard GTF format with transcript and exon information
- **Classification File**: Tab-separated file with transcript classification data

### Output Files

```
output_directory/
├── hub.txt                              # Main hub configuration
├── genomes.txt                          # Genome assembly mapping
├── {genome}_classification.txt          # Copy of classification file
├── {genome}_sqanti3.bb                  # Main bigBed with all transcripts
├── {genome}_sqanti3_full-splice_match.bb      # Category-specific bigBed
├── {genome}_sqanti3_novel_in_catalog.bb       # Category-specific bigBed
├── {genome}_sqanti3_*.bb                      # Other category bigBeds
├── {genome}_sqanti3_*.html                    # Track description pages
├── {output}_{genome}_SQANTI3_Hub.html   # Hub description page
├── {genome}/
│   ├── trackDb.txt                      # Track configuration
│   ├── groups.txt                       # Track groups
│   ├── trix.ix                          # Trix search index
│   └── trix.ixx                         # Trix search index
└── table_reports/                       # (if --tables flag used)
    ├── full-splice_match_isoforms.html
    ├── novel_in_catalog_isoforms.html
    └── ...                              # One HTML per category
```

### Trix Index Files

The `trix.ix` and `trix.ixx` files enable free-text search in the UCSC Genome Browser:
- **trix.ix**: Main index file containing searchable terms
- **trix.ixx**: Secondary index for efficient lookups

These files are automatically generated and linked in `trackDb.txt`.

### Public Hub Readiness Checklist

- Edit `hub.txt` shortLabel/longLabel to be descriptive.
- Ensure `email` is valid and monitored.
- Provide meaningful HTML pages:
  - Top-level hub description (about/usage): `{output}_{genome}_SQANTI3_Hub.html`.
  - Per-track HTML entries (auto-generated; edit as needed for clarity and references).
- Verify `genomes.txt` references `trackDb.txt` and `groups.txt` correctly.
- Host all files at stable URLs (e.g., GitHub Pages or institutional web).
- Test the hub URL in UCSC (My Data → Track Hubs → My Hubs) and confirm filters/search.

## License

All source code is under GNU public license 3.0 (see https://www.gnu.org/licenses/gpl-3.0.de.html).

## Citation

If you use this tool in your research, please cite:

- SQANTI3: https://github.com/ConesaLab/SQANTI3
- UCSC Genome Browser: https://genome.ucsc.edu/
