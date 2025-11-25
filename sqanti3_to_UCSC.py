#!/usr/bin/env python3
"""
SQANTI3 to UCSC Genome Browser Hub Converter

This script converts SQANTI3 output files (*_corrected.gtf and *_classification.txt)
to bigBed format for visualization in the UCSC Genome Browser with hub functionality.

Usage:
    python sqanti3_to_UCSC.py --gtf <gtf_file> --classification <classification_file> --output <output_dir> --genome <genome> [--github-repo username/repository]

Author: Carolina Monzo
"""

import argparse
import os
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path
import pandas as pd
import logging
from collections import defaultdict
# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class SQANTI3ToBigBed:
    def __init__(self, gtf_file, classification_file, output_dir, genome, chrom_sizes_file=None, github_repo=None, github_branch="main", star_sj=None, two_bit_file=None, validate_only=False, dry_run=False):
        self.gtf_file = gtf_file
        self.classification_file = classification_file
        self.output_dir = Path(output_dir)
        self.genome = genome
        self.chrom_sizes_file = chrom_sizes_file
        self.github_repo = github_repo
        self.github_branch = github_branch
        self.temp_dir = None
        self.star_sj = star_sj
        self.star_bigbed = None
        self.two_bit_file = two_bit_file
        self.validate_only = validate_only
        self.dry_run = dry_run
        self.keep_temp = False
        self.category_bigbeds = {}
        
        # Create output directory (handles both relative and absolute paths)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Output directory: {self.output_dir.absolute()}")
        
        # Check if required tools are available
        self._check_dependencies()
    
    def _check_dependencies(self):
        """Check if required UCSC tools are available"""
        required_tools = ['gtfToGenePred', 'genePredToBed', 'bedToBigBed']
        # twoBitInfo required only if using --twobit
        if self.two_bit_file:
            required_tools.append('twoBitInfo')
        missing_tools = []
        
        for tool in required_tools:
            if shutil.which(tool) is None:
                missing_tools.append(tool)
        
        if missing_tools:
            logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            logger.error("Please install UCSC tools running the following script: bash install_ucsc_tools.sh")
            logger.error("Or install them manually from: http://hgdownload.soe.ucsc.edu/admin/exe/")
            logger.error("Or use conda: conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed ucsc-bedtobigbed")
            sys.exit(1)
    
    def _create_temp_dir(self):
        """Create temporary directory for intermediate files"""
        self.temp_dir = tempfile.mkdtemp(prefix="sqanti3_bigbed_")
        logger.info(f"Created temporary directory: {self.temp_dir}")
    
    def _cleanup_temp_dir(self):
        """Clean up temporary directory"""
        if self.temp_dir and os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
            logger.info("Cleaned up temporary directory")
    
    def extract_chrom_sizes(self):
        """Extract chromosome sizes from GTF file"""
        logger.info("Extracting chromosome sizes from GTF file...")
        
        chrom_max_pos = defaultdict(int)
        
        try:
            with open(self.gtf_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if line.startswith('#') or not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 5:
                        continue
                    
                    try:
                        chrom = parts[0]
                        end = int(parts[4])
                        chrom_max_pos[chrom] = max(chrom_max_pos[chrom], end)
                    except (ValueError, IndexError):
                        continue
            
            # Write chrom.sizes file
            chrom_sizes_file = os.path.join(self.temp_dir, "chrom.sizes")
            with open(chrom_sizes_file, 'w') as f:
                for chrom in sorted(chrom_max_pos.keys()):
                    f.write(f"{chrom}\t{chrom_max_pos[chrom]}\n")
            
            logger.info(f"Extracted chromosome sizes for {len(chrom_max_pos)} chromosomes")
            return chrom_sizes_file
            
        except Exception as e:
            logger.error(f"Error extracting chromosome sizes: {e}")
            return None

    def extract_chrom_sizes_from_twobit(self):
        """Use twoBitInfo to compute chromosome sizes from a .2bit file"""
        if not self.two_bit_file:
            return None
        logger.info("Extracting chromosome sizes from 2bit file using twoBitInfo...")
        try:
            # Run twoBitInfo and capture output
            result = subprocess.run(['twoBitInfo', self.two_bit_file, 'stdout'], capture_output=True, text=True, check=True)
            lines = result.stdout.strip().split('\n')
            sizes = []
            for line in lines:
                parts = line.strip().split()
                if len(parts) >= 2:
                    chrom = parts[0]
                    try:
                        size = int(parts[1])
                    except ValueError:
                        continue
                    sizes.append((chrom, size))
            # Sort by size descending (as in UCSC examples). Sorting is not strictly required by bedToBigBed.
            sizes.sort(key=lambda x: x[1], reverse=True)
            chrom_sizes_file = os.path.join(self.temp_dir, 'chrom.sizes')
            with open(chrom_sizes_file, 'w') as out:
                for chrom, size in sizes:
                    out.write(f"{chrom}\t{size}\n")
            logger.info(f"Chrom sizes written from 2bit: {chrom_sizes_file}")
            return chrom_sizes_file
        except subprocess.CalledProcessError as e:
            logger.error(f"twoBitInfo failed: {e}")
            if e.stderr:
                logger.error(e.stderr)
            return None
    
    def parse_classification_file(self):
        """Parse the SQANTI3 classification file and extract filter values"""
        logger.info("Parsing classification file...")
        
        try:
            # Load ALL columns to provide full information
            self.classification_df = pd.read_csv(
                self.classification_file, 
                sep='	',
                dtype={'isoform': 'string'}
            )
            
            # Rename 'isoform' to 'name' for merging
            self.classification_df = self.classification_df.rename(columns={'isoform': 'name'})
            
            logger.info(f"Loaded classification data for {len(self.classification_df)} transcripts with {len(self.classification_df.columns)} columns")
            return True
            
        except Exception as e:
            logger.error(f"Error parsing classification file: {e}")
            return False

    def convert_gtf_to_genepred(self):
        """Convert GTF to GenePred format"""
        logger.info("Converting GTF to GenePred format...")
        
        genepred_file = os.path.join(self.temp_dir, "transcripts.genepred")
        cmd = [
            'gtfToGenePred',
            '-genePredExt',
            '-allErrors',
            '-ignoreGroupsWithoutExons',
            self.gtf_file,
            genepred_file
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info("GTF to GenePred conversion completed")
            return genepred_file
        except subprocess.CalledProcessError as e:
            logger.error(f"GTF to GenePred conversion failed: {e}")
            if e.stderr:
                logger.error(f"Error details: {e.stderr}")
            return None
    
    def convert_genepred_to_bed(self, genepred_file):
        """Convert GenePred to BED format and fix malformed block arrays"""
        logger.info("Converting GenePred to BED format...")
        
        try:
            # First convert with genePredToBed
            temp_bed_file = os.path.join(self.temp_dir, "transcripts_temp.bed")
            cmd = ['genePredToBed', genepred_file, temp_bed_file]
            
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            # Fix the malformed BED file (remove trailing commas from blockSizes and chromStarts)
            bed_file = os.path.join(self.temp_dir, "transcripts.bed")
            with open(temp_bed_file, 'r') as infile, open(bed_file, 'w') as outfile:
                for line in infile:
                    parts = line.strip().split('\t')
                    if len(parts) >= 12:  # BED12 format
                        # Fix blockSizes (field 11) - remove trailing comma
                        if parts[10].endswith(','):
                            parts[10] = parts[10].rstrip(',')
                        
                        # Fix chromStarts (field 12) - remove trailing comma
                        if parts[11].endswith(','):
                            parts[11] = parts[11].rstrip(',')
                        
                        outfile.write('\t'.join(parts) + '\n')
                    else:
                        outfile.write(line)
            
            logger.info("GenePred to BED conversion completed")
            return bed_file
            
        except subprocess.CalledProcessError as e:
            logger.error(f"GenePred to BED conversion failed: {e}")
            if e.stderr:
                logger.error(f"Error details: {e.stderr}")
            return None
    

    
    def _generate_autosql_schema(self, extra_cols, output_path):
        """Generate AutoSQL schema for BED12 + extra columns"""
        with open(output_path, 'w') as f:
            f.write("table sqanti3Transcripts\n")
            f.write('"SQANTI3 transcript annotations with full classification data"\n')
            f.write("(\n")
            f.write('    string chrom;        "Reference sequence chromosome or scaffold"\n')
            f.write('    uint   chromStart;   "Start position in chromosome"\n')
            f.write('    uint   chromEnd;     "End position in chromosome"\n')
            f.write('    string name;         "Transcript ID"\n')
            f.write('    uint   score;        "Score (0-1000)"\n')
            f.write('    char[1] strand;      "+ or -"\n')
            f.write('    uint   thickStart;   "Start of thick display"\n')
            f.write('    uint   thickEnd;     "End of thick display"\n')
            f.write('    uint   itemRgb;      "Item color packed as R*256^2+G*256+B"\n')
            f.write('    int    blockCount;   "Number of blocks (exons)"\n')
            f.write('    int[blockCount] blockSizes;  "Comma separated list of block sizes"\n')
            f.write('    int[blockCount] chromStarts; "Start positions relative to chromStart"\n')
            
            for col in extra_cols:
                # Sanitize column name for AutoSQL (alphanumeric + underscore)
                safe_col = col.replace('.', '_').replace(' ', '_').replace('-', '_').replace('/', '_')
                # Use lstring for safety (handles long text)
                f.write(f'    lstring {safe_col}; "{col}"\n')
            
            f.write(")\n")

    def add_classification_data_to_bed(self, bed_file):
        """Add classification data to BED file using Pandas"""
        logger.info("Adding classification data to BED file (Vectorized)...")
        
        output_file = os.path.join(self.temp_dir, "transcripts_full.bed")
        
        try:
            # Load BED file (headerless)
            bed_cols = [
                'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 
                'blockSizes', 'chromStarts'
            ]
            
            bed_df = pd.read_csv(
                bed_file, 
                sep='	', 
                names=bed_cols,
                dtype={'chrom': 'string', 'name': 'string', 'strand': 'category'}
            )
            
            if not hasattr(self, 'classification_df'):
                logger.error("Classification data not loaded")
                return None

            # Identify extra columns: those in classification but NOT in standard BED12
            # We exclude standard BED columns to avoid duplicates and merge conflicts
            exclude_cols = set(bed_cols)
            exclude_cols.add('ORF_seq')  # User requested to exclude this long column
            
            self.extra_cols = [
                c for c in self.classification_df.columns 
                if c not in exclude_cols and c != 'name'
            ]
            
            # Prepare classification subset for merging
            # We need 'name' for the merge, plus the extra columns
            class_subset = self.classification_df[['name'] + self.extra_cols].copy()
            
            # Merge
            merged = bed_df.merge(class_subset, on='name', how='left')
            
            # Handle categorical data and NaNs
            categorical_defaults = {
                'structural_category': 'NA',
                'subcategory': 'NA',
                'coding': 'NA',
                'FSM_class': 'NA'
            }
            
            for col, default_val in categorical_defaults.items():
                if col in merged.columns:
                    # If it's a categorical column, we must add the category first
                    if isinstance(merged[col].dtype, pd.CategoricalDtype):
                        if default_val not in merged[col].cat.categories:
                            merged[col] = merged[col].cat.add_categories([default_val])
                    
                    merged[col] = merged[col].fillna(default_val)

            # Fill all other extra columns with 'NA' if they have missing values
            for col in self.extra_cols:
                if col in merged.columns:
                     merged[col] = merged[col].fillna('NA')

            # Calculate itemRgb
            def pack_rgb(r, g, b):
                return (r << 16) + (g << 8) + b

            cat_palette = {
                "full-splice_match": pack_rgb(107, 174, 214),
                "incomplete-splice_match": pack_rgb(252, 141, 89),
                "novel_in_catalog": pack_rgb(120, 198, 121),
                "novel_not_in_catalog": pack_rgb(238, 106, 80),
                "genic": pack_rgb(150, 150, 150),
                "antisense": pack_rgb(102, 194, 164),
                "fusion": pack_rgb(218, 165, 32),
                "intergenic": pack_rgb(233, 150, 122),
                "genic_intron": pack_rgb(65, 182, 196),
                "NA": pack_rgb(200, 200, 200)
            }
            
            default_color = pack_rgb(200, 200, 200)
            
            if 'structural_category' in merged.columns:
                merged['itemRgb'] = merged['structural_category'].map(cat_palette).fillna(default_color).astype(int)
            else:
                merged['itemRgb'] = default_color
            
            # Format FSM class
            if 'FSM_class' in merged.columns:
                merged['FSM_class'] = merged['FSM_class'].astype(str)
                mask = merged['FSM_class'] != 'NA'
                merged.loc[mask, 'FSM_class'] = 'FSM' + merged.loc[mask, 'FSM_class']
            
            # Write to file
            final_cols = bed_cols + self.extra_cols
            
            merged[final_cols].to_csv(
                output_file, 
                sep='	', 
                header=False, 
                index=False,
                quoting=3
            )
            
            logger.info(f"Created enhanced BED file with {len(self.extra_cols)} extra columns: {output_file}")
            return output_file
            
        except Exception as e:
            logger.error(f"Error in vectorized BED processing: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None

    def _generate_trix_index(self, enhanced_bed_file, genome_dir):
        """Generate Trix (.ix/.ixx) text index"""
        ixixx_path = shutil.which('ixIxx')
        if not ixixx_path:
            logger.warning("ixIxx not found; skipping Trix index generation")
            return False
            
        logger.info("Generating Trix index...")
        try:
            trix_input = os.path.join(self.temp_dir, 'trix_input.txt')
            
            # We need a set of transcript IDs present in the BED file
            valid_ids = set()
            with open(enhanced_bed_file, 'r') as f:
                for line in f:
                    parts = line.split('\t')
                    if len(parts) > 3:
                        valid_ids.add(parts[3])
            
            # Now stream the classification file and extract keywords for valid IDs
            with open(self.classification_file, 'r') as f_in, open(trix_input, 'w') as f_out:
                header = f_in.readline().rstrip('\n').split('\t')
                try:
                    idx_isoform = header.index('isoform')
                except ValueError:
                    logger.error("isoform column not found in classification file")
                    return False
                    
                skip_cols = ['isoform', 'chrom', 'strand'] 
                useful_cols = [c for c in header if c not in skip_cols]
                col_indices = {c: header.index(c) for c in useful_cols}
                
                for line in f_in:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) <= idx_isoform:
                        continue
                        
                    isoform = parts[idx_isoform]
                    if isoform in valid_ids:
                        keywords = [isoform]
                        for col, idx in col_indices.items():
                            if idx < len(parts):
                                val = parts[idx]
                                if val and val not in ['NA', 'nan', '']:
                                    keywords.append(f"{col}:{val}")
                                    keywords.append(val)
                        
                        unique_keywords = sorted(list(set(keywords)))
                        f_out.write(' '.join(unique_keywords) + '\n')
            
            ix_file = genome_dir / "trix.ix"
            ixx_file = genome_dir / "trix.ixx"
            
            cmd = [ixixx_path, trix_input, str(ix_file), str(ixx_file)]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"Trix index created: {ix_file}, {ixx_file}")
            return True
            
        except Exception as e:
            logger.error(f"Error generating Trix index: {e}")
            return False

    def _create_star_sj_bigbed(self, star_sj_file):
        """Convert STAR splice junctions to bigBed"""
        logger.info("Converting STAR splice junctions to bigBed...")
        try:
            if self.chrom_sizes_file:
                chrom_sizes = self.chrom_sizes_file
            elif self.two_bit_file:
                chrom_sizes = self.extract_chrom_sizes_from_twobit()
            else:
                chrom_sizes = self.extract_chrom_sizes()

            if not chrom_sizes:
                 return None

            bed_file = os.path.join(self.temp_dir, "star_junctions.bed")
            with open(star_sj_file, 'r') as infile, open(bed_file, 'w') as outfile:
                for i, line in enumerate(infile):
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue
                    chrom = parts[0]
                    start = int(parts[1]) - 1
                    end = int(parts[2])
                    strand_val = parts[3]
                    strand = '+' if strand_val == '1' else ('-' if strand_val == '2' else '.')
                    unique = int(parts[6])
                    multi = int(parts[7])
                    score = min(unique + multi, 1000)
                    name = f"JUNC{i+1}"
                    outfile.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
            
            sorted_bed = os.path.join(self.temp_dir, "star_junctions.sorted.bed")
            env = os.environ.copy()
            env["LC_COLLATE"] = "C"
            subprocess.run(['sort', '-k1,1', '-k2,2n', bed_file, '-o', sorted_bed], check=True, env=env)
            
            bb_file = self.output_dir / f"{self.genome}_star_sj.bb"
            subprocess.run(['bedToBigBed', sorted_bed, chrom_sizes, str(bb_file)], check=True)
            logger.info(f"STAR bigBed created: {bb_file}")
            return bb_file
        except Exception as e:
            logger.error(f"Error converting STAR SJ: {e}")
            return None

    def _get_github_raw_url(self, filename, subdir=None):
        """Generate a raw GitHub URL for a file."""
        if self.github_repo:
            base_url = f"https://raw.githubusercontent.com/{self.github_repo}/{self.github_branch}"
            if subdir:
                return f"{base_url}/{subdir}/{filename}"
            else:
                return f"{base_url}/{filename}"
        else:
            # Fallback to relative paths if no GitHub repo is provided
            if subdir:
                return f"{subdir}/{filename}"
            else:
                return filename
    
    def _get_rgb_color(self, structural_category):
        """Get RGB color for structural category - using original hex colors from color_bed.py"""
        # Original hex colors from Alejandro Paniagua's color_bed.py
        cat_palette = {
            "full-splice_match": "6BAED6",
            "incomplete-splice_match": "FC8D59",
            "novel_in_catalog": "78C679",
            "novel_not_in_catalog": "EE6A50",
            "genic": "969696",
            "antisense": "66C2A4",
            "fusion": "DAA520",
            "intergenic": "E9967A",
            "genic_intron": "41B6C4"
        }
        
        # Get HEX code
        hex_color = cat_palette.get(structural_category.lower(), "FFFFFF")
        
        # Convert HEX to RGB (following color_bed.py logic)
        rgb = []
        for i in (0, 2, 4):
            decimal = int(hex_color[i:i+2], 16)
            rgb.append(str(decimal))
        
        return ",".join(rgb)
    


    
    def create_bigbed_file(self, bed_file):
        """Convert full BED to bigBed using dynamic autoSql"""
        logger.info("Converting BED to bigBed...")
        
        try:
            # Determine chromosome sizes
            if self.chrom_sizes_file and os.path.exists(self.chrom_sizes_file):
                chrom_sizes_file = self.chrom_sizes_file
                logger.info(f"Using provided chrom.sizes file: {chrom_sizes_file}")
            elif self.two_bit_file:
                chrom_sizes_file = self.extract_chrom_sizes_from_twobit()
                if not chrom_sizes_file:
                    raise Exception("twoBitInfo failed to generate chrom.sizes")
            else:
                chrom_sizes_file = self.extract_chrom_sizes()
                if not chrom_sizes_file:
                    raise Exception("Could not determine chromosome sizes")

            # bed_file is now transcripts_full.bed (unsorted)
            sorted_bed_file = os.path.join(self.temp_dir, "transcripts_full.sorted.bed")
            env = os.environ.copy()
            env["LC_COLLATE"] = "C"
            sort_cmd = ['sort', '-k1,1', '-k2,2n', bed_file, '-o', sorted_bed_file]
            subprocess.run(sort_cmd, check=True, capture_output=True, text=True, env=env)
            logger.info(f"Sorted BED written to: {sorted_bed_file}")

            bigbed_file = self.output_dir / f"{self.genome}_sqanti3.bb"
            
            # Write dynamic autoSql schema file
            as_path = os.path.join(self.temp_dir, 'sqanti3_schema.as')
            self._generate_autosql_schema(self.extra_cols, as_path)

            # Run bedToBigBed
            # type=bed12+N where N is number of extra columns
            # -tab is required because fields may contain spaces
            num_extra = len(self.extra_cols)
            cmd = ['bedToBigBed', '-tab', f'-as={as_path}', f'-type=bed12+{num_extra}', '-extraIndex=name', sorted_bed_file, chrom_sizes_file, str(bigbed_file)]
            
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"BigBed file created: {bigbed_file}")
            return bigbed_file

        except Exception as e:
            logger.error(f"Error creating bigBed file: {e}")
            if 'result' in locals() and hasattr(result, 'stderr'):
                 logger.error(f"bedToBigBed stderr: {result.stderr}")
            return None

    def create_category_bigbeds(self, main_bed_file):
        """Create separate bigBed files for each structural category"""
        logger.info("Creating category-specific bigBed files...")
        
        try:
            as_path = os.path.join(self.temp_dir, 'sqanti3_schema.as')
            num_extra = len(self.extra_cols)
            chrom_sizes_file = self.chrom_sizes_file or os.path.join(self.temp_dir, 'chrom.sizes')
            
            # Read full BED file to filter by category
            # We use the same columns as written in add_classification_data_to_bed
            bed_cols = [
                'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 
                'blockSizes', 'chromStarts'
            ] + self.extra_cols
            
            df = pd.read_csv(main_bed_file, sep='	', names=bed_cols, dtype={'chrom': 'string'})
            
            if 'structural_category' not in df.columns:
                logger.warning("structural_category column missing, cannot create category tracks")
                return False
                
            categories = df['structural_category'].unique()
            self.category_bigbeds = {}
            
            for cat in categories:
                if pd.isna(cat) or cat == 'NA': 
                    continue
                
                cat_str = str(cat)
                safe_cat_file = cat_str.replace(' ', '_').replace('/', '_')
                
                sub_bed = os.path.join(self.temp_dir, f"sqanti3_{safe_cat_file}.bed")
                sub_sorted = os.path.join(self.temp_dir, f"sqanti3_{safe_cat_file}.sorted.bed")
                sub_bb = self.output_dir / f"{self.genome}_sqanti3_{safe_cat_file}.bb"
                
                # Filter
                cat_df = df[df['structural_category'] == cat]
                
                # Write
                cat_df.to_csv(sub_bed, sep='	', header=False, index=False, quoting=3)
                
                # Sort
                env = os.environ.copy()
                env["LC_COLLATE"] = "C"
                subprocess.run(['sort', '-k1,1', '-k2,2n', sub_bed, '-o', sub_sorted], check=True, env=env)
                
                # bedToBigBed - use -tab because fields may contain spaces
                cmd = ['bedToBigBed', '-tab', f'-as={as_path}', f'-type=bed12+{num_extra}', '-extraIndex=name', sub_sorted, chrom_sizes_file, str(sub_bb)]
                
                result = subprocess.run(cmd, check=False, capture_output=True, text=True)
                if result.returncode != 0:
                    logger.error(f"bedToBigBed failed for category {cat}")
                    logger.error(f"Command: {' '.join(cmd)}")
                    logger.error(f"Stderr: {result.stderr}")
                    raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
                
                self.category_bigbeds[cat_str] = sub_bb.name
                logger.info(f"Created category track: {sub_bb}")
                
            return True
            
        except Exception as e:
            logger.error(f"Error creating category bigBeds: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False

    def create_hub_files(self, bigbed_file):
        """Create UCSC Genome Browser hub files"""
        logger.info("Creating hub files...")
        
        # Extract the output directory name for hub naming
        output_dir_name = self.output_dir.name
        hub_name = f"{output_dir_name}_{self.genome}_SQANTI3_Hub"
        
        # Create hub.txt
        hub_file = self.output_dir / "hub.txt"
        with open(hub_file, 'w', newline='\n') as f:
            f.write(f"hub {hub_name}\n")
            f.write(f"shortLabel {hub_name}\n")
            f.write(f"longLabel SQANTI3 Transcriptome Analysis for {self.genome}\n")
            f.write(f"genomesFile {self._get_github_raw_url('genomes.txt')}\n")
            if self.github_repo:
                f.write(f"email {self.github_repo.split('/')[0]}@users.noreply.github.com\n")
            else:
                f.write(f"email sqanti3_user@users.noreply.github.com\n")
            f.write(f"descriptionUrl {self._get_github_raw_url('README.md')}\n")
        
        # Create genomes.txt
        genomes_file = self.output_dir / "genomes.txt"
        with open(genomes_file, 'w', newline='\n') as f:
            f.write(f"genome {self.genome}\n")
            f.write(f"trackDb {self._get_github_raw_url('trackDb.txt', subdir=self.genome)}\n")
            f.write(f"groups {self._get_github_raw_url('groups.txt', subdir=self.genome)}\n")
        
        # Create genome-specific directory and trackDb.txt
        genome_dir = self.output_dir / self.genome
        genome_dir.mkdir(exist_ok=True)
        
        # Create groups.txt
        groups_file = genome_dir / "groups.txt"
        with open(groups_file, 'w', newline='\n') as gf:
            gf.write("name transcripts\n")
            gf.write("label Transcripts\n")
            gf.write("priority 1\n")
            gf.write("defaultIsClosed 0\n\n")
            gf.write("name junctions\n")
            gf.write("label Splice Junctions\n")
            gf.write("priority 2\n")
            gf.write("defaultIsClosed 0\n")

        # Create per-track HTML files
        transcripts_html_name = f"{self.genome}_sqanti3_track.html"
        transcripts_html = self.output_dir / transcripts_html_name
        with open(transcripts_html, 'w') as tf:
            tf.write(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{self.genome} SQANTI3 Transcripts</title>
    <style>
        body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }}
        h1 {{ color: #333; }}
        ul {{ list-style-type: disc; margin-left: 20px; }}
        code {{ background-color: #f4f4f4; padding: 2px 4px; border-radius: 4px; }}
    </style>
</head>
<body>
    <h1>SQANTI3 Transcripts ({self.genome})</h1>
    <p>This track displays SQANTI3 transcript models. Colors indicate structural category.</p>
    <h3>Filtering</h3>
    <ul>
        <li>Search by name-encoded tokens (e.g., <code>intergenic</code>, <code>mono-exon</code>, <code>FSM1</code>).</li>
        <li>Use per-field filters and ranges (length, exons, coverage, expression) available in the track settings.</li>
    </ul>
    <h3>Color Legend</h3>
    <ul>
        <li>Full-splice match: #6BAED6</li>
        <li>Incomplete-splice match: #FC8D59</li>
        <li>Novel in catalog: #78C679</li>
        <li>Novel not in catalog: #EE6A50</li>
        <li>Genic: #969696</li>
        <li>Antisense: #66C2A4</li>
        <li>Fusion: #DAA520</li>
        <li>Intergenic: #E9967A</li>
        <li>Genic intron: #41B6C4</li>
    </ul>
</body>
</html>""")

        star_html_name = None
        if self.star_bigbed and os.path.exists(self.star_bigbed):
            star_html_name = f"{self.genome}_star_sj_track.html"
            star_html = self.output_dir / star_html_name
            with open(star_html, 'w') as sh:
                sh.write(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{self.genome} STAR Splice Junctions</title>
    <style>
        body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }}
        h1 {{ color: #333; }}
    </style>
</head>
<body>
    <h1>STAR Splice Junctions ({self.genome})</h1>
    <p>Junctions converted from STAR <code>SJ.out.tab</code> (bed6/bigBed).</p>
    <p>Strand and score reflect STAR output; junction name is sequential.</p>
    <p>Sort and range navigation recommended for large regions.</p>
</body>
</html>""")
        
        # Determine number of extra fields
        num_extra = len(self.extra_cols) if hasattr(self, 'extra_cols') else 0

        trackdb_file = genome_dir / "trackDb.txt"
        with open(trackdb_file, 'w', newline='\n') as f:
            f.write(f"track {self.genome}_sqanti3\n")
            f.write(f"bigDataUrl {self._get_github_raw_url(f'{self.genome}_sqanti3.bb')}\n")
            f.write(f"shortLabel SQANTI3 Transcripts\n")
            f.write(f"longLabel SQANTI3 Transcriptome Analysis Results\n")
            f.write(f"type bigBed 12 + {num_extra}\n")
            
            # Add standard filters if columns exist
            def sanitize(col):
                return col.replace('.', '_').replace(' ', '_').replace('-', '_').replace('/', '_')
            
            existing_fields = set([sanitize(c) for c in self.extra_cols])
            
            if 'structural_category' in existing_fields: f.write("filter.structural_category on\n")
            if 'subcategory' in existing_fields: f.write("filter.subcategory on\n")
            if 'coding' in existing_fields: f.write("filter.coding on\n")
            if 'FSM_class' in existing_fields: f.write("filter.FSM_class on\n")
            
            # Range filters
            if 'length' in existing_fields: f.write("filterByRange.length 0:50000\n")
            if 'exons' in existing_fields: f.write("filterByRange.exons 0:100\n")
            if 'min_cov' in existing_fields: f.write("filterByRange.min_cov 0:1000\n")
            if 'iso_exp' in existing_fields: f.write("filterByRange.iso_exp 0:1000\n")
            
            f.write(f"visibility dense\n")
            f.write(f"group transcripts\n")
            f.write(f"itemRgb on\n")
            f.write(f"color 107,174,214\n")
            f.write(f"priority 1\n")
            f.write(f"html {self._get_github_raw_url(transcripts_html_name)}\n")
            f.write(f"searchIndex name\n")
            f.write(f"filterByRange.blockCount 0:100\n")
            f.write(f"filterLabel.blockCount Number of exons\n")
            
            # Trix
            trix_ix = genome_dir / 'trix.ix'
            if os.path.exists(trix_ix):
                f.write(f"searchTrix {self._get_github_raw_url('trix.ix', subdir=self.genome)}\n")

            # STAR
            if self.star_bigbed and os.path.exists(self.star_bigbed):
                f.write("\n")
                f.write(f"track {self.genome}_star_sj\n")
                f.write(f"bigDataUrl {self._get_github_raw_url(self.star_bigbed.name)}\n")
                f.write(f"shortLabel STAR Junctions\n")
                f.write(f"longLabel STAR splice junctions (SJ.out.tab)\n")
                f.write(f"type bigBed 6\n")
                f.write(f"visibility dense\n")
                f.write(f"group junctions\n")
                f.write(f"priority 2\n")
                if star_html_name:
                    f.write(f"html {self._get_github_raw_url(star_html_name)}\n")

            # Category tracks
            if hasattr(self, 'category_bigbeds') and self.category_bigbeds:
                for cat, bb_filename in self.category_bigbeds.items():
                    safe_cat = cat.replace(' ', '_').replace('/', '_')
                    cat_html_name = f"{self.genome}_sqanti3_{safe_cat}.html"
                    
                    # Create HTML
                    script_dir = Path(__file__).parent
                    template_path = script_dir / 'track_template.html'
                    if not template_path.exists():
                        template_path = Path('track_template.html')
                    
                    if template_path.exists():
                        try:
                            with open(template_path, 'r') as f_template:
                                html_content = f_template.read()

                            if hasattr(self, 'classification_df'):
                                count = len(self.classification_df[self.classification_df['structural_category'] == cat])
                            else:
                                count = 0

                            html_content = html_content.replace('%%TITLE%%', f"{self.genome} SQANTI3 {cat} Transcripts")
                            html_content = html_content.replace('%%CATEGORY_NAME%%', cat)
                            html_content = html_content.replace('%%TRANSCRIPT_COUNT%%', str(count))
                            html_content = html_content.replace('%%GENOME%%', self.genome)

                            with open(self.output_dir / cat_html_name, 'w') as ch:
                                ch.write(html_content)
                        except Exception as e:
                            logger.error(f"Error creating HTML for category {cat}: {e}")
                    
                    f.write("\n")
                    track_name = f"{self.genome}_sqanti3_{safe_cat}"
                    short = f"SQANTI3 {cat}"
                    f.write(f"track {track_name}\n")
                    f.write(f"bigDataUrl {self._get_github_raw_url(bb_filename)}\n")
                    f.write(f"shortLabel {short}\n")
                    f.write(f"longLabel SQANTI3 {cat} transcripts\n")
                    f.write(f"type bigBed 12 + {num_extra}\n")
                    f.write(f"visibility hide\n")
                    f.write(f"group transcripts\n")
                    f.write(f"itemRgb on\n")
                    f.write(f"priority 3\n")
                    f.write(f"html {self._get_github_raw_url(cat_html_name)}\n")
                    
                    # Add same filters
                    if 'subcategory' in existing_fields: f.write("filter.subcategory on\n")
                    if 'coding' in existing_fields: f.write("filter.coding on\n")
                    if 'FSM_class' in existing_fields: f.write("filter.FSM_class on\n")
                    if 'length' in existing_fields: f.write("filterByRange.length 0:50000\n")
                    if 'exons' in existing_fields: f.write("filterByRange.exons 0:100\n")
                    if 'min_cov' in existing_fields: f.write("filterByRange.min_cov 0:1000\n")
                    if 'iso_exp' in existing_fields: f.write("filterByRange.iso_exp 0:1000\n")
                    f.write(f"filterByRange.blockCount 0:100\n")
                    f.write(f"filterLabel.blockCount Number of exons\n")

        # README creation (same as before)
        readme_file = self.output_dir / "README.md"
        with open(readme_file, 'w') as f_md:
            f_md.write(f"""# {hub_name}

This hub displays SQANTI3 transcriptome analysis results for the {self.genome} genome assembly.

## 🔍 Advanced Filtering

**This hub uses the bigBed 12+{num_extra} format with native UCSC filters.**

You can filter transcripts by:
- **Structural Category:** FSM, ISM, NIC, NNC, genic, antisense, fusion, intergenic, genic_intron
- **Subcategory:** mono-exon, multi-exon, novel, known, canonical, non-canonical
- **Coding Status:** coding, non_coding, partial_coding, pseudo
- **FSM Class:** A, B, C, D
- **Length:** Range-based filtering for transcript length
- **Coverage:** Range-based filtering for minimum coverage
- **Expression:** Range-based filtering for isoform expression

*Right-click on the track and select "Configure" or "Filter" to access these controls.*

## 🎨 Color Legend

- **Full-splice Match (FSM):** #6BAED6 (Blue)
- **Incomplete-splice Match (ISM):** #FC8D59 (Orange)
- **Novel In Catalog (NIC):** #78C679 (Green)
- **Novel Not In Catalog (NNC):** #EE6A50 (Red)
- **Genic:** #969696 (Gray)
- **Antisense:** #66C2A4 (Teal)
- **Fusion:** #DAA520 (Gold)
- **Intergenic:** #E9967A (Salmon)
- **Genic Intron:** #41B6C4 (Cyan)

## 📊 Classification Data

This hub visualizes data from the SQANTI3 classification analysis.

## 🚀 Usage Instructions

To use this hub in the UCSC Genome Browser:

1. Upload all files in this directory to a web-accessible location (e.g., GitHub).
2. In the UCSC Genome Browser, go to **My Data → Track Hubs**.
3. Enter the URL to your `hub.txt` file.
4. Select the appropriate genome assembly ({self.genome}).
5. The SQANTI3 tracks will appear in your track list.

## 📧 Contact Information

**Note:** The `hub.txt` file includes a placeholder or GitHub-specific email address. If you need to use a different email:
- Edit the `hub.txt` file.
- Change the `email` line to your preferred email address.
- Re-upload the updated file.
""")
        
        logger.info("Hub files created successfully")
        return True
    def run(self):
        """Run the complete conversion pipeline"""
        try:
            self._create_temp_dir()
            
            if not self.parse_classification_file():
                return False

            # If only validating inputs, also check tools and files then exit
            if self.validate_only:
                logger.info("Validation successful: tools present, inputs readable, classification parsed.")
                return True
            
            genepred_file = self.convert_gtf_to_genepred()
            if not genepred_file:
                return False
            
            # Convert GenePred to BED using UCSC tools
            bed_file = self.convert_genepred_to_bed(genepred_file)
            if not bed_file:
                raise Exception("GenePred to BED conversion failed")
            
            # Add classification data to BED file
            bed_file = self.add_classification_data_to_bed(bed_file)
            
            # Exit early after preparing enhanced BED if dry-run requested
            if self.dry_run:
                logger.info("Dry run complete: generated intermediate BED with classification/colors.")
                logger.info(f"Intermediate BED: {bed_file}")
                return True

            # Generate Trix index (after name encoding) if ixIxx is available
            genome_dir = self.output_dir / self.genome
            genome_dir.mkdir(exist_ok=True)
            self._generate_trix_index(bed_file, genome_dir)

            # Create bigBed file
            bigbed_file = self.create_bigbed_file(bed_file)
            if not bigbed_file:
                return False
            
            # Create category-specific tracks
            self.create_category_bigbeds(bed_file)
            
            # Optionally create STAR junctions track
            if self.star_sj:
                self.star_bigbed = self._create_star_sj_bigbed(self.star_sj)
            
            if not self.create_hub_files(bigbed_file):
                return False
            
            logger.info("Conversion completed successfully!")
            logger.info(f"Output directory: {self.output_dir}")
            logger.info(f"BigBed file: {bigbed_file}")
            logger.info("Hub files created and ready for upload to UCSC Genome Browser")
            logger.info("Use searchIndex and filterValues on the name field for comprehensive filtering in the Genome Browser")
            
            return True
            
        except Exception as e:
            logger.error(f"Conversion failed: {e}")
            # Preserve temp directory for debugging
            if self.temp_dir and os.path.exists(self.temp_dir):
                logger.info(f"Temporary directory preserved for debugging: {self.temp_dir}")
            return False
        finally:
            # Cleanup temp directory unless requested to keep it
            if not self.keep_temp:
                self._cleanup_temp_dir()
            else:
                logger.info(f"Temporary directory preserved: {self.temp_dir}")

def main():
    parser = argparse.ArgumentParser(description='Convert SQANTI3 output to UCSC Genome Browser hub for visualization')
    parser.add_argument('--gtf', required=True, help='SQANTI3 corrected GTF file')
    parser.add_argument('--classification', required=True, help='SQANTI3 classification file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--genome', required=True, help='Genome assembly name (e.g., hg38, mm10)')
    parser.add_argument('--chrom-sizes', help='Optional: Path to chromosome sizes file')
    parser.add_argument('--twobit', help='Optional: Genome .2bit file to compute chrom.sizes via twoBitInfo')
    parser.add_argument('--github-repo', help='GitHub repository (format: username/repository) for raw URLs')
    parser.add_argument('--github-branch', default='main', help='GitHub branch (default: main)')
    parser.add_argument('--star-sj', help='Optional: STAR SJ.out.tab to convert into a splice junction track')
    parser.add_argument('--validate-only', action='store_true', help='Validate tools and inputs only, then exit')
    parser.add_argument('--dry-run', action='store_true', help='Prepare intermediates (BED with classification) and exit before bigBed/hub generation')
    parser.add_argument('--keep-temp', action='store_true', help='Keep temporary files for debugging')
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not os.path.exists(args.gtf):
        logger.error(f"GTF file not found: {args.gtf}")
        sys.exit(1)
    
    if not os.path.exists(args.classification):
        logger.error(f"Classification file not found: {args.classification}")
        sys.exit(1)
    
    # Run conversion
    converter = SQANTI3ToBigBed(
        args.gtf,
        args.classification,
        args.output,
        args.genome,
        args.chrom_sizes,
        args.github_repo,
        args.github_branch,
        star_sj=args.star_sj,
        two_bit_file=args.twobit,
        validate_only=args.validate_only,
        dry_run=args.dry_run
    )
    converter.keep_temp = args.keep_temp
    success = converter.run()
    
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
