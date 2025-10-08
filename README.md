# SQANTI-browser

SQANTI-browser converts SQANTI3 transcriptome analysis outputs into UCSC Genome Browser track hubs. It takes SQANTI3 `*_corrected.gtf` and `*_classification.txt`, produces bigBed files, and generates all hub configuration needed for immediate visualization in the UCSC Genome Browser.

## Features

- **Automatic Conversion**: Converts SQANTI3 GTF and classification files to bigBed format
- **Color Coding**: Transcripts are automatically colored by structural category using the same color scheme as SQANTI3
- **Search-Friendly Names**: Encodes classification info in the BED name to enable search/filter in UCSC
- **Hub Generation**: Creates complete UCSC Genome Browser hub files
- **Chromosome Sizes**: Uses provided `--chrom-sizes` file, or extracts sizes directly from the GTF if not provided
- **Comprehensive Logging**: Detailed progress and error reporting

## Prerequisites

### Required Software

1. **UCSC Tools** (required for file format conversion):
   - `gtfToGenePred`
   - `genePredToBed` 
   - `bedToBigBed`

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
   
   # Make executable
   chmod +x gtfToGenePred genePredToBed bedToBigBed
   
   # Add to PATH or move to a directory in your PATH
   sudo mv gtfToGenePred genePredToBed bedToBigBed /usr/local/bin/
   ```
   
   **Option C: Conda installation**
   ```bash
   conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed ucsc-bedtobigbed
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

- **Trix Text Index** (faster text search)
  - Generate a Trix text index by adding `--enable-trix` (requires UCSC `ixIxx` in PATH):
  ```bash
  python sqanti3_to_UCSC.py \
      --gtf your_corrected.gtf \
      --classification your_classification.txt \
      --output output_directory \
      --genome hg38 \
      --enable-trix
  ```
  - This produces `genome/trix.ix` and `genome/trix.ixx` and wires `searchTrix trix.ix` in `trackDb.txt`.

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

- **bigBed 12+ with fielded filters**
  - Enable fielded filters and numeric ranges with `--bb12plus`:
  ```bash
  python sqanti3_to_UCSC.py \
      --gtf your_corrected.gtf \
      --classification your_classification.txt \
      --output output_directory \
      --genome hg38 \
      --bb12plus
  ```
  - Adds per-field `filter` and `filterByRange` controls for structural codes, subcategory, coding, FSM, length, exons, coverage, expression.

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
| `--enable-trix` | No | Generate Trix (.ix/.ixx) text index if `ixIxx` is available |
| `--star-sj` | No | Path to STAR `SJ.out.tab` to add a splice junctions track |

**Output Directory Flexibility:**
- **Relative paths**: `./my_project`, `../analysis`
- **Absolute paths**: `~/Documents/research`, `/home/username/projects`
- **Home directory**: `~/my_hub`, `$HOME/transcriptome_analysis`

**GitHub Integration:**
- **Automatic raw URLs**: When `--github-repo` is specified, all file references use GitHub raw URLs via `https://raw.githubusercontent.com/{repo}/{branch}/...`
- **No manual editing**: Hub files are ready to use directly in UCSC Genome Browser
- **Branch support**: Specify different branches with `--github-branch`

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

## Filter Values

All columns from your SQANTI3 classification file are preserved as filter values, allowing you to:

- Filter transcripts by structural category
- Filter by associated gene
- Filter by any other classification criteria
- Use the UCSC Genome Browser's built-in filtering interface

## Uploading to UCSC Genome Browser

### Step 1: Host Your Files

Upload all generated files to a web-accessible location (e.g., GitHub Pages, your institution's web server, or cloud storage with public URLs).

### Step 2: Add the Hub

1. Go to the [UCSC Genome Browser](https://genome.ucsc.edu/)
2. Navigate to **My Data** → **Track Hubs**
3. Enter the URL to your `hub.txt` file
4. Click **Add Hub**

### Step 3: Select Genome Assembly

1. Choose the appropriate genome assembly from the dropdown
2. Your SQANTI3 tracks will appear in the track list
3. Enable the tracks you want to view

## Example Workflow

### 1. Run SQANTI3

```bash
# Your existing SQANTI3 pipeline
python SQANTI3.py -g reference.gtf -j reference.gff3 -o output_dir input.fasta
```

### 2. Run the converter

```bash
python sqanti3_to_UCSC.py \
    --gtf output_dir/input_corrected.gtf \
    --classification output_dir/input_classification.txt \
    --output hub_output \
    --genome hg38
```

### 3. Upload and View

1. Upload the `hub_output` directory to your web server (for example, your GitHub repository)
2. Check that the paths in your `hub.txt` file correspond to where your files are accessible in the web
3. Go to https://genome.ucsc.edu/index.html
4. Click on My Data -> Track Hubs -> My Hubs and add the link to your publicly available hub.txt file into the URL window.
4. Click on Add Hub.
5. Navigate to your region of interest
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

- **bigBed**: Binary format for efficient display of large datasets
   - `{genome}_sqanti3.bb` - The main bigBed file containing transcript data
- **Hub Files**: UCSC Genome Browser configuration files
   - `hub.txt` - Main hub configuration file
   - `genomes.txt` - Genome assembly mapping
   - `{genome}/trackDb.txt` - Track configuration for the specific genome
   - `{output_directory}_{genome}_SQANTI3_Hub.html` - Detailed description and usage instructions

**Note**: The hub name includes the output directory name for uniqueness (e.g., `my_project_hg38_SQANTI3_Hub`)

```
output_directory/
├── hg38_sqanti3.bb
├── hub.txt
├── genomes.txt
├── hg38/
│   └── trackDb.txt
└── hg38_SQANTI3_Hub.html
```

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