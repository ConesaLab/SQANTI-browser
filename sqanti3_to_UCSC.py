#!/usr/bin/env python3
"""
SQANTI3 to UCSC Genome Browser Hub Converter

This script converts SQANTI3 output files (*_corrected.gtf and *_classification.txt)
to bigBed format for visualization in the UCSC Genome Browser with hub functionality.

Usage:
    python sqanti3_to_UCSC.py --gtf <gtf_file> --classification <classification_file> --output <output_dir> --genome <genome>

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

# Try to import the report generator
try:
    from filter_isoforms import generate_html_reports
except ImportError:
    generate_html_reports = None

class SQANTI3ToBigBed:
    # Valid sort-by options for isoform ordering
    VALID_SORT_OPTIONS = [
        'iso_exp',      # Default: isoform expression (highest first)
        'length',       # Transcript length (longest first)
        'FL',           # Full-length reads (highest first)
        'diff_to_TSS',  # Distance to reference TSS
        'diff_to_TTS',  # Distance to reference TTS
        'diff_to_gene_TSS',   # Distance to gene TSS
        'diff_to_gene_TTS',   # Distance to gene TTS
        'dist_to_CAGE_peak',  # Distance to CAGE peak
        'dist_to_polyA_site'  # Distance to polyA site
    ]
    
    def __init__(self, gtf_file, classification_file, output_dir, genome, chrom_sizes_file=None, star_sj=None, cage_peaks=None, polya_peaks=None, ref_gtf=None, two_bit_file=None, validate_only=False, dry_run=False, sort_by='iso_exp', no_category_tracks=False):
        self.gtf_file = gtf_file
        self.classification_file = classification_file
        self.output_dir = Path(output_dir)
        self.genome = genome
        self.chrom_sizes_file = chrom_sizes_file
        self.temp_dir = None
        self.star_sj = star_sj
        self.cage_peaks = cage_peaks
        self.polya_peaks = polya_peaks
        self.ref_gtf = ref_gtf
        self.star_bigbed = None
        self.cage_bigbed = None
        self.polya_bigbed = None
        self.ref_bigbed = None
        self.no_category_tracks = no_category_tracks
        self.two_bit_file = two_bit_file
        self.validate_only = validate_only
        self.dry_run = dry_run
        self.keep_temp = False
        self.category_bigbeds = {}
        self.sort_by = sort_by
        
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

            # Calculate itemRgb as packed integer (R*65536 + G*256 + B)
            # This matches the 'uint itemRgb' definition in the .as schema
            def pack_rgb(r, g, b):
                return r * 65536 + g * 256 + b
            
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

    def _sort_bed_for_visualization(self, df, sort_by=None):
        """Sort BED DataFrame for optimal visualization in UCSC Genome Browser.
        
        Sorting strategy:
        1. Primary: chromosome and start position (REQUIRED by bedToBigBed format)
        2. Secondary: associated_transcript (groups isoforms at same position together)
        3. Tertiary: sort_by metric (default: iso_exp, highest first)
        
        Note: bigBed format requires strict (chrom, chromStart) sorting. Isoforms of
        the same reference transcript but with different start positions will appear
        at their respective genomic positions. Within the same start position, isoforms
        are grouped by reference transcript and sorted by the chosen metric.
        
        Args:
            df: pandas DataFrame with BED columns
            sort_by: Column to sort by within each reference transcript group.
                     If None, uses self.sort_by. Options: iso_exp, length, FL,
                     diff_to_TSS, diff_to_TTS, diff_to_gene_TSS, diff_to_gene_TTS,
                     dist_to_CAGE_peak, dist_to_polyA_site
        
        Returns:
            Sorted DataFrame
        """
        if sort_by is None:
            sort_by = self.sort_by
        
        if sort_by == 'none':
            sort_by = None
        
        # Check if required columns exist
        has_associated_transcript = 'associated_transcript' in df.columns
        
        if not has_associated_transcript:
            logger.warning("Column 'associated_transcript' not found. Isoforms will not be grouped by reference transcript.")
        
        # Determine the sort metric column
        original_sort_by = sort_by
        fallback_used = False
        
        # Check if sort column exists and has valid (non-NA) values
        if sort_by is not None:
            if sort_by not in df.columns:
                logger.warning(f"Sort column '{sort_by}' not found in data.")
                sort_by = None
            elif df[sort_by].isna().all() or (df[sort_by].astype(str) == 'NA').all():
                logger.warning(f"Column '{sort_by}' contains only NA values.")
                sort_by = None
            else:
                # Check for partial NA values
                na_count = df[sort_by].isna().sum() + (df[sort_by].astype(str) == 'NA').sum()
                if na_count > 0:
                    logger.warning(f"Column '{sort_by}' has {na_count} NA values out of {len(df)} rows. NA values will be sorted last.")
        
        # Fallback logic - only if user requested a specific sort that failed
        if sort_by is None and original_sort_by is not None and original_sort_by != 'none':
            # We don't try to fallback to iso_exp anymore, we just fallback to no sorting
            logger.warning(f"Could not sort by '{original_sort_by}'. Falling back to genomic position order.")

        # Create a numeric version of the sort column for proper sorting
        if sort_by:
            sort_col_numeric = f'_sort_key_{sort_by}'
            df[sort_col_numeric] = pd.to_numeric(df[sort_by], errors='coerce')
        
        # Build sort parameters
        # bigBed REQUIRES sorting by (chrom, chromStart) - this cannot be changed
        # Within the same position, we can group by transcript and sort by metric
        sort_cols = ['chrom', 'chromStart']
        ascending = [True, True]
        
        if has_associated_transcript:
            sort_cols.append('associated_transcript')
            ascending.append(True)
        
        if sort_by:
            sort_cols.append(sort_col_numeric)
            # For distance metrics, smaller (closer) might be better, so ascending
            # For expression/length/FL, higher is typically better, so descending
            if sort_by in ['diff_to_TSS', 'diff_to_TTS', 'diff_to_gene_TSS', 'diff_to_gene_TTS', 
                           'dist_to_CAGE_peak', 'dist_to_polyA_site']:
                # These are distances - sort ascending (closest first)
                ascending.append(True)
            else:
                # iso_exp, length, FL - higher is better, sort descending
                ascending.append(False)
        
        # Perform the sort
        df_sorted = df.sort_values(
            by=sort_cols,
            ascending=ascending,
            na_position='last',
            key=lambda col: col if col.name != 'chrom' else col.astype(str)
        )
        
        # Clean up temporary sort column
        if sort_by:
            df_sorted = df_sorted.drop(columns=[sort_col_numeric])
        
        sort_info = f"by {sort_by}" if sort_by else "by genomic position only"
        group_info = "grouped by reference transcript, " if has_associated_transcript else ""
        logger.info(f"Sorted {len(df_sorted)} transcripts: {group_info}{sort_info}")
        
        return df_sorted

    def _generate_trix_index(self, enhanced_bed_file, genome_dir):
        """Generate Trix (.ix/.ixx) text index with rich descriptions for fast search
        
        This matches the working format from October 30 version:
        - TAB-separated with 3 columns: ID, description, synonyms
        - Uses classification_df for field values
        """
        ixixx_path = shutil.which('ixIxx')
        if not ixixx_path:
            logger.warning("ixIxx not found; skipping Trix index generation")
            return False

        logger.info("Generating Trix index...")
        try:
            trix_input = os.path.join(self.temp_dir, 'trix_input.txt')
            
            # Build a dictionary from classification_df for fast lookup
            if not hasattr(self, 'classification_df'):
                logger.error("Classification data not loaded")
                return False
            
            classification_data = self.classification_df.set_index('name').to_dict('index')
            
            with open(enhanced_bed_file, 'r') as bed_in, open(trix_input, 'w') as t_out:
                for line in bed_in:
                    parts = line.rstrip('\t\n').split('\t')
                    if len(parts) >= 12:
                        tid = parts[3]
                        d = classification_data.get(tid, {})
                        
                        # Build field dictionary
                        fields = {
                            'structural_category': d.get('structural_category', 'unknown'),
                            'chrom': d.get('chrom', ''),
                            'strand': d.get('strand', ''),
                            'subcategory': d.get('subcategory', 'unknown'),
                            'coding': d.get('coding', 'unknown'),
                            'FSM_class': d.get('FSM_class', '0'),
                            'length': d.get('length', '0'),
                            'exons': d.get('exons', '0'),
                            'min_cov': d.get('min_cov', '0'),
                            'iso_exp': d.get('iso_exp', '0'),
                            'associated_gene': d.get('associated_gene', ''),
                            'associated_transcript': d.get('associated_transcript', ''),
                            'ref_length': d.get('ref_length', ''),
                            'ref_exons': d.get('ref_exons', ''),
                            'diff_to_TSS': d.get('diff_to_TSS', ''),
                            'diff_to_TTS': d.get('diff_to_TTS', ''),
                            'diff_to_gene_TSS': d.get('diff_to_gene_TSS', ''),
                            'diff_to_gene_TTS': d.get('diff_to_gene_TTS', ''),
                            'RTS_stage': d.get('RTS_stage', ''),
                            'all_canonical': d.get('all_canonical', ''),
                            'min_sample_cov': d.get('min_sample_cov', ''),
                            'min_cov_pos': d.get('min_cov_pos', ''),
                            'sd_cov': d.get('sd_cov', ''),
                            'FL': d.get('FL', ''),
                            'n_indels': d.get('n_indels', ''),
                            'n_indels_junc': d.get('n_indels_junc', ''),
                            'bite': d.get('bite', ''),
                            'gene_exp': d.get('gene_exp', ''),
                            'ratio_exp': d.get('ratio_exp', ''),
                            'ORF_length': d.get('ORF_length', ''),
                            'CDS_length': d.get('CDS_length', ''),
                            'CDS_start': d.get('CDS_start', ''),
                            'CDS_end': d.get('CDS_end', ''),
                            'CDS_genomic_start': d.get('CDS_genomic_start', ''),
                            'CDS_genomic_end': d.get('CDS_genomic_end', ''),
                            'predicted_NMD': d.get('predicted_NMD', ''),
                            'perc_A_downstream_TTS': d.get('perc_A_downstream_TTS', ''),
                            'seq_A_downstream_TTS': d.get('seq_A_downstream_TTS', ''),
                            'dist_to_CAGE_peak': d.get('dist_to_CAGE_peak', ''),
                            'within_CAGE_peak': d.get('within_CAGE_peak', ''),
                            'dist_to_polyA_site': d.get('dist_to_polyA_site', ''),
                            'within_polyA_site': d.get('within_polyA_site', ''),
                            'polyA_motif': d.get('polyA_motif', ''),
                            'polyA_dist': d.get('polyA_dist', ''),
                            'polyA_motif_found': d.get('polyA_motif_found', ''),
                            'ratio_TSS': d.get('ratio_TSS', '')
                        }
                        
                        # Build description: "tid key value key value ..."
                        desc = f"{tid} " + ' '.join([f"{k} {v}" for k, v in fields.items() if v not in (None, '', 'NA', 'nan')])
                        
                        # Build synonyms: values + prefixed versions for all columns
                        # This allows searching like "strand_plus" or "structural_category_intergenic"
                        synonym_parts = []
                        for k, v in fields.items():
                            if v not in (None, '', 'NA', 'nan'):
                                v_str = str(v)
                                synonym_parts.append(v_str)
                                # For prefixed version, convert standalone + and - (strand values)
                                # ixIxx strips these characters when they're standalone
                                if v_str == '+':
                                    v_safe = 'plus'
                                elif v_str == '-':
                                    v_safe = 'minus'
                                else:
                                    v_safe = v_str
                                # Add prefixed version (e.g., strand_plus, structural_category_incomplete-splice_match)
                                # Note: colon (:) doesn't work in UCSC search, use underscore
                                synonym_parts.append(f"{k}_{v_safe}")
                        synonyms = ' '.join(synonym_parts)
                        
                        # Write TAB-separated: ID \t description \t synonyms
                        t_out.write(f"{tid}\t{desc}\t{synonyms}\n")

            ix_path = os.path.join(genome_dir, 'trix.ix')
            ixx_path = os.path.join(genome_dir, 'trix.ixx')
            # Use -maxWordLength=64 to handle long prefixed terms like
            # structural_category_novel_not_in_catalog (40 chars, default is 31)
            cmd = ['ixIxx', '-maxWordLength=64', trix_input, ix_path, ixx_path]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"Trix index generated: {ix_path}, {ixx_path}")
            return True
            
        except Exception as e:
            logger.warning(f"Failed to generate Trix index: {e}")
            import traceback
            logger.warning(traceback.format_exc())
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
            
            bb_file = self.output_dir / self.genome / f"{self.genome}_star_sj.bb"
            subprocess.run(['bedToBigBed', sorted_bed, chrom_sizes, str(bb_file)], check=True)
            logger.info(f"STAR bigBed created: {bb_file}")
            return bb_file
        except Exception as e:
            logger.error(f"Error converting STAR SJ: {e}")
            return None

    def _create_cage_bigbed(self, cage_file):
        """Convert CAGE peaks BED to bigBed"""
        logger.info("Converting CAGE peaks to bigBed...")
        try:
            if self.chrom_sizes_file:
                chrom_sizes = self.chrom_sizes_file
            elif self.two_bit_file:
                chrom_sizes = self.extract_chrom_sizes_from_twobit()
            else:
                chrom_sizes = self.extract_chrom_sizes()

            if not chrom_sizes:
                 return None

            # Input BED format: chr, start, end, name, score, strand, ...
            # We want to create a standard BED6 or BED9+ file
            # The example file has 9 columns: chr, start, end, name, score, strand, thickStart, thickEnd, itemRgb
            
            # Create a sorted BED file
            bed_file = os.path.join(self.temp_dir, "cage_peaks.bed")
            
            # Read valid chromosomes and their sizes from chrom.sizes
            valid_chroms = {}
            with open(chrom_sizes, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        valid_chroms[parts[0]] = int(parts[1])
            
            # Filter and sort the input file
            # We assume the input is already in a valid BED format (BED3+)
            # We sort by chrom, start
            env = os.environ.copy()
            env["LC_COLLATE"] = "C"
            
            # Filter cage_file to only include chromosomes in valid_chroms AND within bounds
            filtered_bed = os.path.join(self.temp_dir, "cage_peaks_filtered.bed")
            filtered_count = 0
            total_count = 0
            skipped_bounds = 0
            
            with open(cage_file, 'r') as infile, open(filtered_bed, 'w') as outfile:
                for line in infile:
                    total_count += 1
                    parts = line.split('\t')
                    if parts and len(parts) >= 3 and parts[0] in valid_chroms:
                        try:
                            # Check bounds: end position must be <= chrom size
                            # BED end is 0-based exclusive, chrom size is 1-based length
                            # So end coordinate can be at most chrom_size
                            chrom = parts[0]
                            end_pos = int(parts[2])
                            
                            if end_pos <= valid_chroms[chrom]:
                                outfile.write(line)
                                filtered_count += 1
                            else:
                                skipped_bounds += 1
                        except ValueError:
                            # Skip lines with non-integer coordinates
                            continue
            
            if filtered_count == 0:
                logger.warning(f"No CAGE peaks kept after filtering. Total: {total_count}. Skipped out of bounds: {skipped_bounds}")
                return None
            
            if filtered_count < total_count:
                logger.warning(f"Filtered CAGE peaks to match genome: {filtered_count}/{total_count} peaks kept. ({skipped_bounds} out of bounds, {total_count - filtered_count - skipped_bounds} on unknown chromosomes)")

            # Read first line to determine number of columns to set bigBed type
            num_cols = 0
            with open(filtered_bed, 'r') as f:
                line = f.readline()
                if line:
                    num_cols = len(line.strip().split('\t'))
            
            if num_cols < 3:
                logger.error("CAGE file must have at least 3 columns")
                return None
                
            sorted_bed = os.path.join(self.temp_dir, "cage_peaks.sorted.bed")
            subprocess.run(['sort', '-k1,1', '-k2,2n', filtered_bed, '-o', sorted_bed], check=True, env=env)
            
            bb_file = self.output_dir / self.genome / f"{self.genome}_cage_peaks.bb"
            
            # Determine type argument (e.g., bed9)
            # Default to bed6 if uncertain, but example suggests bed9
            # bedToBigBed will fail if type doesn't match column count
            type_arg = f"bed{min(num_cols, 12)}"
            
            # Special handling for bed9: columns 7-9 are thickStart, thickEnd, itemRgb
            # If the file has 9 columns, bedToBigBed expects them to be thickStart, thickEnd, itemRgb
            if num_cols == 9:
                type_arg = "bed9"
            elif num_cols > 9:
                # If > 9 cols, we specify bed9+N
                type_arg = f"bed9+{num_cols-9}"
            
            # If using AutoSQL, we could be more specific, but for now we rely on standard BED types
            
            # For the specific example file: human.refTSS_v3.1.hg38.bed
            # It has 9 columns. UCSC bedToBigBed needs correct type.
            # If it has 9 cols, we use bed9
            
            subprocess.run(['bedToBigBed', f'-type={type_arg}', '-tab', sorted_bed, chrom_sizes, str(bb_file)], check=True)
            logger.info(f"CAGE bigBed created: {bb_file}")
            return bb_file
        except Exception as e:
            logger.error(f"Error converting CAGE peaks: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None

    def _create_polya_bigbed(self, polya_file):
        """Convert PolyA peaks BED to bigBed"""
        logger.info("Converting PolyA peaks to bigBed...")
        try:
            if self.chrom_sizes_file:
                chrom_sizes = self.chrom_sizes_file
            elif self.two_bit_file:
                chrom_sizes = self.extract_chrom_sizes_from_twobit()
            else:
                chrom_sizes = self.extract_chrom_sizes()

            if not chrom_sizes:
                 return None
            
            # Create a sorted BED file
            # Read valid chromosomes and their sizes from chrom.sizes
            valid_chroms = {}
            with open(chrom_sizes, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        valid_chroms[parts[0]] = int(parts[1])
            
            # Filter and sort the input file
            # We assume the input is already in a valid BED format (BED3+)
            # We sort by chrom, start
            env = os.environ.copy()
            env["LC_COLLATE"] = "C"
            
            # Filter polya_file to only include chromosomes in valid_chroms AND within bounds
            filtered_bed = os.path.join(self.temp_dir, "polya_peaks_filtered.bed")
            filtered_count = 0
            total_count = 0
            skipped_bounds = 0
            
            with open(polya_file, 'r') as infile, open(filtered_bed, 'w') as outfile:
                for line in infile:
                    total_count += 1
                    parts = line.split('\t')
                    if parts and len(parts) >= 3 and parts[0] in valid_chroms:
                        try:
                            # Check bounds: end position must be <= chrom size
                            # BED end is 0-based exclusive, chrom size is 1-based length
                            # So end coordinate can be at most chrom_size
                            chrom = parts[0]
                            end_pos = int(parts[2])
                            
                            if end_pos <= valid_chroms[chrom]:
                                outfile.write(line)
                                filtered_count += 1
                            else:
                                skipped_bounds += 1
                        except ValueError:
                            # Skip lines with non-integer coordinates
                            continue
            
            if filtered_count == 0:
                logger.warning(f"No PolyA peaks kept after filtering. Total: {total_count}. Skipped out of bounds: {skipped_bounds}")
                return None
            
            if filtered_count < total_count:
                logger.warning(f"Filtered PolyA peaks to match genome: {filtered_count}/{total_count} peaks kept. ({skipped_bounds} out of bounds, {total_count - filtered_count - skipped_bounds} on unknown chromosomes)")

            # PolyA peaks often have float scores in column 5, but BED format requires integers (0-1000)
            # Convert the score column to integers
            corrected_bed = os.path.join(self.temp_dir, "polya_peaks_corrected.bed")
            with open(filtered_bed, 'r') as infile, open(corrected_bed, 'w') as outfile:
                for line in infile:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) >= 5:
                        try:
                            # Convert float score to integer (0-1000 range)
                            # Option 1: multiply by 1000 and cap at 1000
                            score_float = float(parts[4])
                            score_int = min(int(score_float * 1000), 1000)
                            parts[4] = str(score_int)
                        except ValueError:
                            # If conversion fails, set to 0
                            parts[4] = '0'
                    outfile.write('\t'.join(parts) + '\n')

            # Read first line to determine number of columns to set bigBed type
            num_cols = 0
            with open(corrected_bed, 'r') as f:
                line = f.readline()
                if line:
                    num_cols = len(line.strip().split('\t'))
            
            if num_cols < 3:
                logger.error("PolyA file must have at least 3 columns")
                return None
                
            sorted_bed = os.path.join(self.temp_dir, "polya_peaks.sorted.bed")
            subprocess.run(['sort', '-k1,1', '-k2,2n', corrected_bed, '-o', sorted_bed], check=True, env=env)
            
            bb_file = self.output_dir / self.genome / f"{self.genome}_polya_peaks.bb"
            
            # PolyA peaks file format is actually BED6+ with extra fields, not true BED9
            # Columns: chrom, start, end, name, score, strand, [extra fields...]
            # Use bed6+N instead of bed9+N
            if num_cols <= 6:
                type_arg = f"bed{num_cols}"
            else:
                # BED6 + extra fields
                type_arg = f"bed6+{num_cols-6}"
            
            subprocess.run(['bedToBigBed', f'-type={type_arg}', '-tab', sorted_bed, chrom_sizes, str(bb_file)], check=True)
            logger.info(f"PolyA bigBed created: {bb_file}")
            return bb_file
        except Exception as e:
            logger.error(f"Error converting PolyA peaks: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None

    def _create_reference_bigbed(self, ref_gtf_file):
        """Convert reference GTF to bigBed for direct comparison"""
        logger.info("Converting reference GTF to bigBed...")
        try:
            # Get chromosome sizes
            if self.chrom_sizes_file:
                chrom_sizes = self.chrom_sizes_file
            elif self.two_bit_file:
                chrom_sizes = self.extract_chrom_sizes_from_twobit()
            else:
                chrom_sizes = self.extract_chrom_sizes()

            if not chrom_sizes:
                return None

            # Convert GTF to GenePred
            ref_genepred = os.path.join(self.temp_dir, "reference.genepred")
            logger.info("Converting reference GTF to GenePred...")
            subprocess.run(['gtfToGenePred', ref_gtf_file, ref_genepred], check=True)
            
            # Convert GenePred to BED12
            ref_bed = os.path.join(self.temp_dir, "reference.bed")
            logger.info("Converting reference GenePred to BED...")
            subprocess.run(['genePredToBed', ref_genepred, ref_bed], check=True)
            
            # Read valid chromosomes and their sizes from chrom.sizes
            valid_chroms = {}
            with open(chrom_sizes, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        valid_chroms[parts[0]] = int(parts[1])
            
            # Filter BED to only include transcripts within valid chromosome bounds
            filtered_bed = os.path.join(self.temp_dir, "reference_filtered.bed")
            filtered_count = 0
            total_count = 0
            skipped_bounds = 0
            skipped_chroms = 0
            
            with open(ref_bed, 'r') as infile, open(filtered_bed, 'w') as outfile:
                for line in infile:
                    total_count += 1
                    parts = line.split('\t')
                    if parts and len(parts) >= 3:
                        try:
                            chrom = parts[0]
                            start_pos = int(parts[1])
                            end_pos = int(parts[2])
                            
                            if chrom not in valid_chroms:
                                skipped_chroms += 1
                                continue
                            
                            # Check if both start and end are within bounds
                            if end_pos <= valid_chroms[chrom] and start_pos >= 0:
                                outfile.write(line)
                                filtered_count += 1
                            else:
                                skipped_bounds += 1
                        except (ValueError, IndexError):
                            continue
            
            if filtered_count == 0:
                logger.warning(f"No reference transcripts kept after filtering. Total: {total_count}. Skipped: {skipped_bounds} out of bounds, {skipped_chroms} on unknown chromosomes")
                return None
            
            if filtered_count < total_count:
                logger.warning(f"Filtered reference transcripts to match genome: {filtered_count}/{total_count} transcripts kept. ({skipped_bounds} out of bounds, {skipped_chroms} on unknown chromosomes)")
            
            # Sort the BED file
            env = os.environ.copy()
            env["LC_COLLATE"] = "C"
            ref_sorted_bed = os.path.join(self.temp_dir, "reference.sorted.bed")
            subprocess.run(['sort', '-k1,1', '-k2,2n', filtered_bed, '-o', ref_sorted_bed], check=True, env=env)
            
            # Convert to bigBed
            bb_file = self.output_dir / self.genome / f"{self.genome}_reference.bb"
            subprocess.run(['bedToBigBed', '-type=bed12', '-tab', ref_sorted_bed, chrom_sizes, str(bb_file)], check=True)
            
            logger.info(f"Reference bigBed created: {bb_file}")
            return bb_file
        except Exception as e:
            logger.error(f"Error converting reference GTF: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None

    def _get_relative_path(self, filename, subdir=None):
        """Generate a relative path for a file in the hub structure."""
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
    
    def _get_category_hex_color(self, structural_category):
        """Get hex color for structural category"""
        cat_palette = {
            "full-splice_match": "#6BAED6",
            "incomplete-splice_match": "#FC8D59",
            "novel_in_catalog": "#78C679",
            "novel_not_in_catalog": "#EE6A50",
            "genic": "#969696",
            "antisense": "#66C2A4",
            "fusion": "#DAA520",
            "intergenic": "#E9967A",
            "genic_intron": "#41B6C4"
        }
        return cat_palette.get(structural_category.lower().replace(' ', '_'), "#6BAED6")
    
    def _generate_filter_options_html(self, category_filter=None):
        """Generate HTML documentation for filter options in classification file column order.
        
        Args:
            category_filter: Optional structural category to filter by (for category-specific tracks)
        
        Returns:
            HTML string with filter options documentation
        """
        if not hasattr(self, 'classification_df'):
            return "<p>No classification data available.</p>"
        
        df = self.classification_df
        if category_filter:
            df = df[df['structural_category'] == category_filter]
        
        # Filters in classification file column order (matching trackDb.txt order)
        # Each entry: ('text', col, label) or ('range', col, min, max, label)
        ordered_filters = [
            ('range', 'length', 0, 100000, 'Transcript length (bp)'),
            ('range', 'exons', 0, 200, 'Number of exons'),
            ('text', 'structural_category', 'Structural Category'),
            ('text', 'associated_gene', 'Associated Gene'),
            ('text', 'associated_transcript', 'Associated Transcript'),
            ('range', 'ref_length', 0, 100000, 'Reference transcript length'),
            ('range', 'ref_exons', 0, 200, 'Reference exon count'),
            ('range', 'diff_to_TSS', -100000, 100000, 'Distance to reference TSS'),
            ('range', 'diff_to_TTS', -100000, 100000, 'Distance to reference TTS'),
            ('range', 'diff_to_gene_TSS', -100000, 100000, 'Distance to gene TSS'),
            ('range', 'diff_to_gene_TTS', -100000, 100000, 'Distance to gene TTS'),
            ('text', 'subcategory', 'Subcategory'),
            ('text', 'RTS_stage', 'RTS Stage'),
            ('text', 'all_canonical', 'All Canonical Splice Sites'),
            ('range', 'min_sample_cov', 0, 10000, 'Minimum sample coverage'),
            ('range', 'min_cov', 0, 10000, 'Minimum junction coverage'),
            ('text', 'min_cov_pos', 'Minimum Coverage Position'),
            ('range', 'sd_cov', 0, 1000, 'Coverage standard deviation'),
            ('range', 'FL', 0, 100000, 'Full-length read count'),
            ('range', 'n_indels', 0, 100, 'Number of indels'),
            ('range', 'n_indels_junc', 0, 100, 'Number of indels at junctions'),
            ('text', 'bite', 'BITE Status'),
            ('range', 'iso_exp', 0, 100000, 'Isoform expression (TPM)'),
            ('range', 'gene_exp', 0, 100000, 'Gene expression (TPM)'),
            ('range', 'ratio_exp', 0, 1, 'Expression ratio (isoform/gene)'),
            ('text', 'FSM_class', 'FSM Class'),
            ('text', 'coding', 'Coding Status'),
            ('range', 'ORF_length', 0, 50000, 'ORF length (aa)'),
            ('range', 'CDS_length', 0, 50000, 'CDS length (bp)'),
            ('range', 'CDS_start', 0, 100000, 'CDS start position'),
            ('range', 'CDS_end', 0, 100000, 'CDS end position'),
            ('text', 'predicted_NMD', 'Predicted NMD'),
            ('range', 'perc_A_downstream_TTS', 0, 100, 'Percent A downstream of TTS'),
            ('range', 'dist_to_CAGE_peak', -10000, 10000, 'Distance to CAGE peak'),
            ('text', 'within_CAGE_peak', 'Within CAGE Peak'),
            ('range', 'dist_to_polyA_site', -10000, 10000, 'Distance to polyA site'),
            ('range', 'polyA_dist', -1000, 1000, 'PolyA distance'),
            ('text', 'polyA_motif_found', 'PolyA Motif Found'),
            ('range', 'ratio_TSS', 0, 10, 'TSS ratio'),
        ]
        
        html_parts = []
        
        # Columns that use dropdown menus
        dropdown_cols = {
            'structural_category', 'subcategory', 'coding', 'FSM_class', 
            'all_canonical', 'RTS_stage', 'bite', 'predicted_NMD',
            'within_CAGE_peak', 'polyA_motif_found'
        }
        
        # Build table with all filters in order
        html_parts.append('<table>')
        html_parts.append('<tr><th>Filter</th><th>Type</th><th>Values / Range</th></tr>')
        
        for filter_def in ordered_filters:
            if filter_def[0] == 'text':
                _, col, label = filter_def
                if col in df.columns:
                    if col in ['associated_gene', 'associated_transcript', 'min_cov_pos']:
                        # These have too many unique values - use text filter
                        unique_count = df[col].nunique()
                        html_parts.append(f'<tr><td><strong>{label}</strong></td><td>🔤 Text</td><td>{unique_count} unique values. Use wildcards like <code>*</code></td></tr>')
                    elif col in dropdown_cols:
                        # Dropdown filter - show values vertically
                        value_counts = df[col].value_counts(dropna=False)
                        unique_vals = []
                        for val, count in value_counts.items():
                            if pd.notna(val) and str(val) != 'NA' and str(val) != '':
                                unique_vals.append(f"<code>{val}</code> ({count})")
                        
                        if unique_vals:
                            # Display vertically with line breaks
                            display_vals = '<br>'.join(unique_vals)
                            html_parts.append(f'<tr><td><strong>{label}</strong></td><td>📋 Dropdown</td><td class="filter-values">{display_vals}</td></tr>')
                    else:
                        # Text filter with values shown
                        value_counts = df[col].value_counts(dropna=False)
                        unique_vals = []
                        for val, count in value_counts.items():
                            if pd.notna(val) and str(val) != 'NA' and str(val) != '':
                                unique_vals.append(f"<code>{val}</code> ({count})")
                        
                        if unique_vals:
                            if len(unique_vals) > 8:
                                display_vals = ', '.join(unique_vals[:8]) + f', ... ({len(unique_vals)} total)'
                            else:
                                display_vals = ', '.join(unique_vals)
                            html_parts.append(f'<tr><td><strong>{label}</strong></td><td>🔤 Text</td><td class="filter-values">{display_vals}</td></tr>')
            else:
                # Range filter
                _, col, min_val, max_val, label = filter_def
                if col in df.columns:
                    # Get actual min/max from data
                    try:
                        data_min = df[col].dropna().astype(float).min()
                        data_max = df[col].dropna().astype(float).max()
                        html_parts.append(f'<tr><td><strong>{label}</strong></td><td>📊 Slider</td><td>Data range: {data_min:.0f} to {data_max:.0f}</td></tr>')
                    except (ValueError, TypeError):
                        html_parts.append(f'<tr><td><strong>{label}</strong></td><td>📊 Slider</td><td>Range: {min_val} to {max_val}</td></tr>')
        
        # Add blockCount at the end
        html_parts.append('<tr><td><strong>Number of exons (from BED)</strong></td><td>📊 Slider</td><td>Range: 0 to 200</td></tr>')
        
        html_parts.append('</table>')
        
        return '\n    '.join(html_parts)

    
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

            # Read, sort with visualization logic, and write
            bed_cols = [
                'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 
                'blockSizes', 'chromStarts'
            ] + self.extra_cols
            
            df = pd.read_csv(bed_file, sep='\t', names=bed_cols, dtype={'chrom': 'string'})
            
            # Apply visualization-optimized sorting
            df_sorted = self._sort_bed_for_visualization(df, self.sort_by)
            
            # Write sorted BED file
            sorted_bed_file = os.path.join(self.temp_dir, "transcripts_full.sorted.bed")
            df_sorted.to_csv(sorted_bed_file, sep='\t', header=False, index=False, quoting=3)
            logger.info(f"Sorted BED written to: {sorted_bed_file}")

            bigbed_file = self.output_dir / self.genome / f"{self.genome}_sqanti3.bb"
            
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
                
                sub_sorted = os.path.join(self.temp_dir, f"sqanti3_{safe_cat_file}.sorted.bed")
                sub_bb = self.output_dir / self.genome / f"{self.genome}_sqanti3_{safe_cat_file}.bb"
                
                # Filter by category
                cat_df = df[df['structural_category'] == cat].copy()
                
                # Apply visualization-optimized sorting
                cat_df_sorted = self._sort_bed_for_visualization(cat_df, self.sort_by)
                
                # Write sorted BED
                cat_df_sorted.to_csv(sub_sorted, sep='\t', header=False, index=False, quoting=3)
                
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
            f.write(f"genomesFile {self._get_relative_path('genomes.txt')}\n")
            f.write(f"email sqanti3_user@users.noreply.github.com\n")
            f.write(f"descriptionUrl {self._get_relative_path('README.md')}\n")
        
        # Create genomes.txt
        genomes_file = self.output_dir / "genomes.txt"
        with open(genomes_file, 'w', newline='\n') as f:
            f.write(f"genome {self.genome}\n")
            f.write(f"trackDb {self._get_relative_path('trackDb.txt', subdir=self.genome)}\n")
            f.write(f"groups {self._get_relative_path('groups.txt', subdir=self.genome)}\n")
            
            # If using 2bit file, add required fields for user-defined genome
            if self.two_bit_file:
                # Copy 2bit file to genome directory
                genome_2bit_name = f"{self.genome}.2bit"
                shutil.copy(self.two_bit_file, self.output_dir / self.genome / genome_2bit_name)
                
                f.write(f"twoBitPath {self._get_relative_path(genome_2bit_name, subdir=self.genome)}\n")
                f.write(f"organism user-defined\n")
                f.write(f"description user-defined\n")
                f.write(f"scientificName user-defined\n")
                f.write(f"orderKey 4800\n")
                
                # Create HTML description for 2bit genome
                genome_html_name = f"{self.genome}_genome.html"
                with open(self.output_dir / self.genome / genome_html_name, 'w') as gf:
                    gf.write(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{self.genome} Description</title>
    <style>
        body {{ font-family: sans-serif; margin: 20px; }}
        h1 {{ color: #006699; }}
        .info-box {{ background: #f0f0f0; padding: 15px; border-radius: 5px; }}
    </style>
</head>
<body>
    <h1>User-Defined Genome: {self.genome}</h1>
    <div class="info-box">
        <p><strong>Organism:</strong> User-defined</p>
        <p><strong>Description:</strong> Custom reference genome provided via 2bit file.</p>
        <p><strong>Source File:</strong> {os.path.basename(self.two_bit_file)}</p>
    </div>
    <h2>Chromosomes</h2>
    <ul>""")
                    # Add chromosome list from chrom.sizes
                    try:
                        with open(os.path.join(self.temp_dir, 'chrom.sizes'), 'r') as cf:
                            for line in cf:
                                parts = line.strip().split('\t')
                                if len(parts) >= 2:
                                    gf.write(f"<li><strong>{parts[0]}:</strong> {parts[1]} bp</li>\n")
                    except Exception:
                        gf.write("<li>Could not read chromosome sizes.</li>")
                    
                    gf.write("""    </ul>
</body>
</html>""")
                
                f.write(f"htmlPath {self._get_relative_path(genome_html_name, subdir=self.genome)}\n")
                
                # Get default position from chrom.sizes (first chromosome)
                try:
                    with open(os.path.join(self.temp_dir, 'chrom.sizes'), 'r') as cf:
                        first_line = cf.readline().strip()
                        if first_line:
                            chrom = first_line.split('\t')[0]
                            f.write(f"defaultPos {chrom}:1-5000\n")
                except Exception as e:
                    logger.warning(f"Could not determine default position for genomes.txt: {e}")
        
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
            gf.write("name validation\n")
            gf.write("label Validation Data\n")
            gf.write("priority 2\n")
            gf.write("defaultIsClosed 0\n")

        # Build filter options documentation from classification data
        filter_options_html = self._generate_filter_options_html()

        # Create per-track HTML files
        transcripts_html_name = f"{self.genome}_sqanti3_track.html"
        transcripts_html = self.output_dir / self.genome / transcripts_html_name
        with open(transcripts_html, 'w') as tf:
            tf.write(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{self.genome} SQANTI3 Transcripts</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif; line-height: 1.6; margin: 20px; background-color: #f8f9fa; }}
        .container {{ max-width: 900px; margin: 0 auto; background: white; padding: 25px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.05); }}
        h1 {{ color: #333; border-bottom: 2px solid #6BAED6; padding-bottom: 10px; }}
        h2 {{ color: #495057; margin-top: 25px; }}
        h3 {{ color: #6c757d; margin-top: 20px; }}
        ul {{ list-style-type: disc; margin-left: 20px; }}
        code {{ background-color: #e9ecef; padding: 2px 6px; border-radius: 4px; font-size: 0.9em; }}
        .color-box {{ display: inline-block; width: 16px; height: 16px; margin-right: 8px; vertical-align: middle; border-radius: 3px; }}
        table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
        th, td {{ border: 1px solid #dee2e6; padding: 10px; text-align: left; }}
        th {{ background-color: #e9ecef; font-weight: 600; }}
        .filter-values {{ font-family: monospace; font-size: 0.85em; color: #495057; }}
        .tip {{ background-color: #d1ecf1; border-left: 4px solid #17a2b8; padding: 10px 15px; margin: 15px 0; border-radius: 0 4px 4px 0; }}
    </style>
</head>
<body>
<div class="container">
    <h1>SQANTI3 Transcripts ({self.genome})</h1>
    <p>This track displays <strong>{len(self.classification_df)}</strong> SQANTI3 transcript models. Colors indicate structural category.</p>
    
    <h2>How to Filter</h2>
    <div class="tip">
        <strong>Tip:</strong> Right-click on the track and select "Configure" to access all filters. 
        Use <code>*</code> as a wildcard in text filters (e.g., <code>FSM*</code> matches FSMA, FSMB, FSMC, FSMD).
    </div>
    
    <h2>Available Filters</h2>
    <p>All filters available for this track:</p>
    
{filter_options_html}
    
    <h2>Color Legend</h2>
    <table>
        <tr><th>Category</th><th>Color</th></tr>
        <tr><td><span class="color-box" style="background-color: #6BAED6;"></span>Full-splice match (FSM)</td><td>#6BAED6</td></tr>
        <tr><td><span class="color-box" style="background-color: #FC8D59;"></span>Incomplete-splice match (ISM)</td><td>#FC8D59</td></tr>
        <tr><td><span class="color-box" style="background-color: #78C679;"></span>Novel in catalog (NIC)</td><td>#78C679</td></tr>
        <tr><td><span class="color-box" style="background-color: #EE6A50;"></span>Novel not in catalog (NNC)</td><td>#EE6A50</td></tr>
        <tr><td><span class="color-box" style="background-color: #969696;"></span>Genic</td><td>#969696</td></tr>
        <tr><td><span class="color-box" style="background-color: #66C2A4;"></span>Antisense</td><td>#66C2A4</td></tr>
        <tr><td><span class="color-box" style="background-color: #DAA520;"></span>Fusion</td><td>#DAA520</td></tr>
        <tr><td><span class="color-box" style="background-color: #E9967A;"></span>Intergenic</td><td>#E9967A</td></tr>
        <tr><td><span class="color-box" style="background-color: #41B6C4;"></span>Genic intron</td><td>#41B6C4</td></tr>
    </table>
</div>
</body>
</html>""")

        star_html_name = None
        if self.star_bigbed and os.path.exists(self.star_bigbed):
            star_html_name = f"{self.genome}_star_sj_track.html"
            star_html = self.output_dir / self.genome / star_html_name
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

        cage_html_name = None
        if self.cage_bigbed and os.path.exists(self.cage_bigbed):
            cage_html_name = f"{self.genome}_cage_peaks_track.html"
            cage_html = self.output_dir / self.genome / cage_html_name
            with open(cage_html, 'w') as ch:
                ch.write(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{self.genome} CAGE Peaks</title>
    <style>
        body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }}
        h1 {{ color: #333; }}
    </style>
</head>
<body>
    <h1>CAGE Peaks ({self.genome})</h1>
    <p>Cap Analysis of Gene Expression (CAGE) peaks, representing transcription start sites (TSS).</p>
    <p>These peaks are used to validate the 5' ends of transcripts.</p>
</body>
</html>""")

        polya_html_name = None
        if self.polya_bigbed and os.path.exists(self.polya_bigbed):
            polya_html_name = f"{self.genome}_polya_peaks_track.html"
            polya_html = self.output_dir / self.genome / polya_html_name
            with open(polya_html, 'w') as ph:
                ph.write(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{self.genome} PolyA Peaks</title>
    <style>
        body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }}
        h1 {{ color: #333; }}
    </style>
</head>
<body>
    <h1>PolyA Peaks ({self.genome})</h1>
    <p>PolyA site peaks, representing transcription termination sites (TTS).</p>
    <p>These peaks are used to validate the 3' ends of transcripts.</p>
</body>
</html>""")

        ref_html_name = None
        if self.ref_bigbed and os.path.exists(self.ref_bigbed):
            ref_html_name = f"{self.genome}_reference_track.html"
            ref_html = self.output_dir / self.genome / ref_html_name
            with open(ref_html, 'w') as rh:
                rh.write(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{self.genome} Reference Annotation</title>
    <style>
        body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }}
        h1 {{ color: #333; }}
    </style>
</head>
<body>
    <h1>Reference Annotation ({self.genome})</h1>
    <p>Reference genome annotation used for SQANTI3 analysis.</p>
    <p>This track displays the reference transcripts for direct comparison with SQANTI3-annotated isoforms.</p>
    <p><strong>Source:</strong> {os.path.basename(self.ref_gtf)}</p>
</body>
</html>""")
        
        # Determine number of extra fields
        num_extra = len(self.extra_cols) if hasattr(self, 'extra_cols') else 0
        
        trackdb_file = genome_dir / "trackDb.txt"
        with open(trackdb_file, 'w', newline='\n') as f:
            f.write(f"track {self.genome}_sqanti3\n")
            # Files are in the same directory as trackDb.txt
            f.write(f"bigDataUrl {self._get_relative_path(f'{self.genome}_sqanti3.bb')}\n")
            f.write(f"shortLabel SQANTI3 Transcripts\n")
            f.write(f"longLabel SQANTI3 Transcriptome Analysis Results\n")
            f.write(f"type bigBed 12 + {num_extra}\n")
            
            # Add comprehensive filters based on LRGASP hub example
            # Filters are ordered to match the classification file column order
            def sanitize(col):
                return col.replace('.', '_').replace(' ', '_').replace('-', '_').replace('/', '_')
            
            existing_fields = set([sanitize(c) for c in self.extra_cols])
            
            # Combined filters in classification file column order
            # Each entry is either:
            #   ('text', col, label) for categorical/text filters
            #   ('range', col, min, max, label) for numeric range filters
            ordered_filters = [
                # First columns from classification file
                ('text', 'strand', 'Strand (+ or -)'),
                ('range', 'length', 0, 100000, 'Transcript length (bp)'),
                ('range', 'exons', 0, 200, 'Number of exons'),
                ('text', 'structural_category', 'Structural category (FSM, ISM, NIC, NNC, etc.)'),
                ('text', 'associated_gene', 'Associated gene'),
                ('text', 'associated_transcript', 'Associated transcript'),
                ('range', 'ref_length', 0, 100000, 'Reference transcript length'),
                ('range', 'ref_exons', 0, 200, 'Reference exon count'),
                ('range', 'diff_to_TSS', -100000, 100000, 'Distance to reference TSS'),
                ('range', 'diff_to_TTS', -100000, 100000, 'Distance to reference TTS'),
                ('range', 'diff_to_gene_TSS', -100000, 100000, 'Distance to gene TSS'),
                ('range', 'diff_to_gene_TTS', -100000, 100000, 'Distance to gene TTS'),
                ('text', 'subcategory', 'Subcategory (mono-exon, multi-exon, etc.)'),
                ('text', 'RTS_stage', 'RTS stage'),
                ('text', 'all_canonical', 'All canonical splice sites'),
                ('range', 'min_sample_cov', 0, 10000, 'Minimum sample coverage'),
                ('range', 'min_cov', 0, 10000, 'Minimum junction coverage'),
                ('text', 'min_cov_pos', 'Minimum coverage position'),
                ('range', 'sd_cov', 0, 1000, 'Coverage standard deviation'),
                ('range', 'FL', 0, 100000, 'Full-length read count'),
                ('range', 'n_indels', 0, 100, 'Number of indels'),
                ('range', 'n_indels_junc', 0, 100, 'Number of indels at junctions'),
                ('text', 'bite', 'BITE status'),
                ('range', 'iso_exp', 0, 100000, 'Isoform expression (TPM)'),
                ('range', 'gene_exp', 0, 100000, 'Gene expression (TPM)'),
                ('range', 'ratio_exp', 0, 1, 'Expression ratio (isoform/gene)'),
                ('text', 'FSM_class', 'FSM class (A, B, C, D)'),
                ('text', 'coding', 'Coding status (coding, non_coding)'),
                ('range', 'ORF_length', 0, 50000, 'ORF length (aa)'),
                ('range', 'CDS_length', 0, 50000, 'CDS length (bp)'),
                ('range', 'CDS_start', 0, 100000, 'CDS start position'),
                ('range', 'CDS_end', 0, 100000, 'CDS end position'),
                ('text', 'predicted_NMD', 'Predicted NMD'),
                ('range', 'perc_A_downstream_TTS', 0, 100, 'Percent A downstream of TTS'),
                ('range', 'dist_to_CAGE_peak', -10000, 10000, 'Distance to CAGE peak'),
                ('text', 'within_CAGE_peak', 'Within CAGE peak'),
                ('range', 'dist_to_polyA_site', -10000, 10000, 'Distance to polyA site'),
                ('range', 'polyA_dist', -1000, 1000, 'PolyA distance'),
                ('text', 'polyA_motif_found', 'PolyA motif found'),
                ('range', 'ratio_TSS', 0, 10, 'TSS ratio'),
            ]
            
            # Columns that should use dropdown (filterValues) - those with few unique values
            # Columns with many unique values (gene, transcript) will use text search
            dropdown_cols = {
                'structural_category', 'subcategory', 'coding', 'FSM_class', 
                'all_canonical', 'RTS_stage', 'bite', 'predicted_NMD',
                'within_CAGE_peak', 'polyA_motif_found'
            }
            
            # Get unique values for dropdown columns from classification data
            def get_unique_values(col):
                if col in self.classification_df.columns:
                    vals = self.classification_df[col].dropna().astype(str).unique()
                    # Filter out NA/empty and sort
                    vals = sorted([v for v in vals if v and v != 'NA' and v != 'nan'])
                    return vals
                return []
            
            # Write filters in order
            for filter_def in ordered_filters:
                if filter_def[0] == 'text':
                    _, col, label = filter_def
                    if col in existing_fields:
                        if col in dropdown_cols:
                            # Use filterValues for dropdown
                            unique_vals = get_unique_values(col)
                            if unique_vals:
                                vals_str = ','.join(unique_vals)
                                f.write(f"filterValues.{col} {vals_str}\n")
                                f.write(f"filterLabel.{col} {label}\n")
                            else:
                                # Fallback to text filter if no values found
                                f.write(f"filterText.{col} *\n")
                                f.write(f"filterType.{col} wildcard\n")
                                f.write(f"filterLabel.{col} {label}\n")
                        else:
                            # Use text filter for columns with many values
                            f.write(f"filterText.{col} *\n")
                            f.write(f"filterType.{col} wildcard\n")
                            f.write(f"filterLabel.{col} {label}\n")
                else:  # range
                    _, col, min_val, max_val, label = filter_def
                    if col in existing_fields:
                        f.write(f"filter.{col} {min_val}:{max_val}\n")
                        f.write(f"filterByRange.{col} on\n")
                        f.write(f"filterLimits.{col} {min_val}:{max_val}\n")
                        f.write(f"filterLabel.{col} {label}\n")
            
            # Block count (exons from BED12) - always available at end
            f.write(f"filter.blockCount 0:200\n")
            f.write(f"filterByRange.blockCount on\n")
            f.write(f"filterLimits.blockCount 0:200\n")
            f.write(f"filterLabel.blockCount Number of exons (from BED)\n")
            
            # Keep references to filter definitions for category tracks
            categorical_filters = {f[1]: f[2] for f in ordered_filters if f[0] == 'text'}
            numeric_filters = {f[1]: (f[2], f[3], f[4]) for f in ordered_filters if f[0] == 'range'}
            
            f.write(f"visibility full\n")
            f.write(f"group transcripts\n")
            f.write(f"itemRgb on\n")
            f.write(f"priority 1\n")
            f.write(f"html {self._get_relative_path(transcripts_html_name)}\n")
            
            # Search Configuration
            f.write(f"searchIndex name\n")
            
            # Add Trix search index if present (relative path, not URL)
            trix_ix = genome_dir / 'trix.ix'
            if os.path.exists(trix_ix):
                f.write(f"searchTrix trix.ix\n")

            # STAR
            if self.star_bigbed and os.path.exists(self.star_bigbed):
                f.write("\n")
                f.write(f"track {self.genome}_star_sj\n")
                # self.star_bigbed is an absolute Path object, we need just the filename
                star_bb_name = os.path.basename(self.star_bigbed)
                f.write(f"bigDataUrl {self._get_relative_path(star_bb_name)}\n")
                f.write(f"shortLabel STAR Junctions\n")
                f.write(f"longLabel STAR splice junctions (SJ.out.tab)\n")
                f.write(f"type bigBed 6\n")
                f.write(f"visibility full\n")
                f.write(f"group validation\n")
                f.write(f"priority 2\n")
                if star_html_name:
                    f.write(f"html {self._get_relative_path(star_html_name)}\n")

            # CAGE
            if self.cage_bigbed and os.path.exists(self.cage_bigbed):
                f.write("\n")
                f.write(f"track {self.genome}_cage_peaks\n")
                cage_bb_name = os.path.basename(self.cage_bigbed)
                f.write(f"bigDataUrl {self._get_relative_path(cage_bb_name)}\n")
                f.write(f"shortLabel CAGE Peaks\n")
                f.write(f"longLabel CAGE Peaks (TSS Validation)\n")
                # Use bed9 if colors are present, otherwise bed6 or whatever was detected
                # The file type isn't strictly checked by hubCheck here but bedToBigBed cares
                # Ideally we pass the type determined during creation, but for now we'll assume standard
                f.write(f"type bigBed\n") 
                f.write(f"visibility dense\n")
                f.write(f"group validation\n")
                f.write(f"priority 3\n")
                f.write(f"color 0,128,0\n") # Green color for CAGE peaks
                if cage_html_name:
                    f.write(f"html {self._get_relative_path(cage_html_name)}\n")

            # PolyA
            if self.polya_bigbed and os.path.exists(self.polya_bigbed):
                f.write("\n")
                f.write(f"track {self.genome}_polya_peaks\n")
                polya_bb_name = os.path.basename(self.polya_bigbed)
                f.write(f"bigDataUrl {self._get_relative_path(polya_bb_name)}\n")
                f.write(f"shortLabel PolyA Peaks\n")
                f.write(f"longLabel PolyA Site Peaks (TTS Validation)\n")
                f.write(f"type bigBed\n") 
                f.write(f"visibility dense\n")
                f.write(f"group validation\n")
                f.write(f"priority 4\n")
                f.write(f"color 200,0,0\n") # Dark red for PolyA peaks
                if polya_html_name:
                    f.write(f"html {self._get_relative_path(polya_html_name)}\n")

            # Reference GTF
            if self.ref_bigbed and os.path.exists(self.ref_bigbed):
                f.write("\n")
                f.write(f"track {self.genome}_reference\n")
                ref_bb_name = os.path.basename(self.ref_bigbed)
                f.write(f"bigDataUrl {self._get_relative_path(ref_bb_name)}\n")
                f.write(f"shortLabel Reference\n")
                f.write(f"longLabel Reference Genome Annotation\n")
                f.write(f"type bigBed 12\n") 
                f.write(f"visibility pack\n")
                f.write(f"group validation\n")
                f.write(f"priority 5\n")
                f.write(f"color 70,70,70\n") # Gray for reference
                if ref_html_name:
                    f.write(f"html {self._get_relative_path(ref_html_name)}\n")

            # Category tracks
            if hasattr(self, 'category_bigbeds') and self.category_bigbeds:
                for cat, bb_filename in self.category_bigbeds.items():
                    safe_cat = cat.replace(' ', '_').replace('/', '_')
                    cat_html_name = f"{self.genome}_sqanti3_{safe_cat}.html"
                    
                    # Get count and filter options for this category
                    if hasattr(self, 'classification_df'):
                        count = len(self.classification_df[self.classification_df['structural_category'] == cat])
                        cat_filter_options = self._generate_filter_options_html(category_filter=cat)
                    else:
                        count = 0
                        cat_filter_options = "<p>No classification data available.</p>"
                    
                    # Create comprehensive HTML for category track
                    cat_html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{self.genome} SQANTI3 {cat} Transcripts</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif; line-height: 1.6; margin: 20px; background-color: #f8f9fa; }}
        .container {{ max-width: 900px; margin: 0 auto; background: white; padding: 25px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.05); }}
        h1 {{ color: #333; border-bottom: 2px solid {self._get_category_hex_color(cat)}; padding-bottom: 10px; }}
        h2 {{ color: #495057; margin-top: 25px; }}
        h3 {{ color: #6c757d; margin-top: 20px; }}
        code {{ background-color: #e9ecef; padding: 2px 6px; border-radius: 4px; font-size: 0.9em; }}
        table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
        th, td {{ border: 1px solid #dee2e6; padding: 10px; text-align: left; }}
        th {{ background-color: #e9ecef; font-weight: 600; }}
        .filter-values {{ font-family: monospace; font-size: 0.85em; color: #495057; }}
        .tip {{ background-color: #d1ecf1; border-left: 4px solid #17a2b8; padding: 10px 15px; margin: 15px 0; border-radius: 0 4px 4px 0; }}
        .category-badge {{ display: inline-block; background-color: {self._get_category_hex_color(cat)}; color: white; padding: 5px 12px; border-radius: 15px; font-weight: bold; }}
    </style>
    </head>
<body>
<div class="container">
    <h1>SQANTI3 Transcripts: <span class="category-badge">{cat}</span></h1>
    <p>This track contains <strong>{count}</strong> transcripts classified as <strong>{cat}</strong> for the <strong>{self.genome}</strong> genome assembly.</p>
    
    <h2>How to Filter</h2>
    <div class="tip">
        <strong>Tip:</strong> Right-click on the track and select "Configure" to access all filters. 
        Use <code>*</code> as a wildcard in text filters (e.g., <code>mono*</code> matches mono-exon).
    </div>
    
    <h2>Available Filters</h2>
    <p>All filters available for this category:</p>
    
    {cat_filter_options}
</div>
</body>
</html>"""
                    
                    try:
                        with open(self.output_dir / self.genome / cat_html_name, 'w') as ch:
                            ch.write(cat_html_content)
                    except Exception as e:
                        logger.error(f"Error creating HTML for category {cat}: {e}")

                    f.write("\n")
                    track_name = f"{self.genome}_sqanti3_{safe_cat}"
                    short = f"SQANTI3 {cat}"
                    f.write(f"track {track_name}\n")
                    # bb_filename is an absolute Path object or string, we need just the filename
                    bb_name = os.path.basename(bb_filename)
                    f.write(f"bigDataUrl {self._get_relative_path(bb_name)}\n")
                    f.write(f"shortLabel {short}\n")
                    f.write(f"longLabel SQANTI3 {cat} transcripts\n")
                    f.write(f"type bigBed 12 + {num_extra}\n")
                    f.write(f"visibility full\n")
                    f.write(f"group transcripts\n")
                    f.write(f"itemRgb on\n")
                    f.write(f"priority 3\n")
                    f.write(f"html {self._get_relative_path(cat_html_name)}\n")
                    
                    # Add comprehensive filters (same order as main track)
                    # Get unique values for this category's data
                    cat_df = self.classification_df[self.classification_df['structural_category'] == cat]
                    
                    def get_cat_unique_values(col):
                        if col in cat_df.columns:
                            vals = cat_df[col].dropna().astype(str).unique()
                            vals = sorted([v for v in vals if v and v != 'NA' and v != 'nan'])
                            return vals
                        return []
                    
                    for filter_def in ordered_filters:
                        if filter_def[0] == 'text':
                            _, col, label = filter_def
                            if col in existing_fields:
                                if col in dropdown_cols:
                                    unique_vals = get_cat_unique_values(col)
                                    if unique_vals:
                                        vals_str = ','.join(unique_vals)
                                        f.write(f"filterValues.{col} {vals_str}\n")
                                        f.write(f"filterLabel.{col} {label}\n")
                                    else:
                                        f.write(f"filterText.{col} *\n")
                                        f.write(f"filterType.{col} wildcard\n")
                                        f.write(f"filterLabel.{col} {label}\n")
                                else:
                                    f.write(f"filterText.{col} *\n")
                                    f.write(f"filterType.{col} wildcard\n")
                                    f.write(f"filterLabel.{col} {label}\n")
                        else:  # range
                            _, col, min_val, max_val, label = filter_def
                            if col in existing_fields:
                                f.write(f"filter.{col} {min_val}:{max_val}\n")
                                f.write(f"filterByRange.{col} on\n")
                                f.write(f"filterLimits.{col} {min_val}:{max_val}\n")
                                f.write(f"filterLabel.{col} {label}\n")
                    
                    f.write(f"filter.blockCount 0:200\n")
                    f.write(f"filterByRange.blockCount on\n")
                    f.write(f"filterLimits.blockCount 0:200\n")
                    f.write(f"filterLabel.blockCount Number of exons (from BED)\n")
        
        # README creation (same as before)
        readme_file = self.output_dir / "README.md"
        with open(readme_file, 'w') as f_md:
            f_md.write(f"""# {hub_name}

This hub displays SQANTI3 transcriptome analysis results for the {self.genome} genome assembly.
This hub visualizes data from the SQANTI3 classification analysis for the {hub_name} sample.

## 🚀 Usage Instructions

To use this hub in the UCSC Genome Browser:

1. Upload all files in this directory to a web-accessible location (e.g., GitHub).
2. In the UCSC Genome Browser, go to **My Data → Track Hubs**.
3. Enter the URL to your `hub.txt` file.
4. Select the appropriate genome assembly ({self.genome}).
5. The SQANTI3 tracks will appear in your track list.

**Important:** The `hub.txt` file includes a placeholder or GitHub-specific email address. If you need to use a different email:
- Edit the `hub.txt` file.
- Change the `email` line to your preferred email address.
- Re-upload the updated file.

## 🔍 Advanced Filtering

**This hub uses the bigBed 12+{num_extra} format with native UCSC filters.**

You can filter transcripts using dropdown menus and range sliders:
- **Structural Category:** full-splice_match, incomplete-splice_match, novel_in_catalog, novel_not_in_catalog, genic, antisense, fusion, intergenic, genic_intron
- **Subcategory:** mono-exon, multi-exon, 3prime_fragment, 5prime_fragment, reference_match, etc.
- **Coding Status:** coding, non_coding
- **FSM Class:** A, B, C, D
- **Boolean Filters:** RTS_stage, bite, predicted_NMD, within_CAGE_peak, polyA_motif_found
- **Numeric Ranges:** length, exons, expression, coverage, ORF_length, and more
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
            
            # Create category-specific tracks (unless disabled)
            if not self.no_category_tracks:
                self.create_category_bigbeds(bed_file)
            else:
                logger.info("Skipping category-specific tracks (--no-category-tracks flag set)")
            
            # Optionally create STAR junctions track
            if self.star_sj:
                self.star_bigbed = self._create_star_sj_bigbed(self.star_sj)
            
            # Optionally create CAGE peaks track
            if self.cage_peaks:
                self.cage_bigbed = self._create_cage_bigbed(self.cage_peaks)
            
            # Optionally create PolyA peaks track
            if self.polya_peaks:
                self.polya_bigbed = self._create_polya_bigbed(self.polya_peaks)
            
            # Optionally create reference GTF track
            if self.ref_gtf:
                self.ref_bigbed = self._create_reference_bigbed(self.ref_gtf)
            
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
    parser.add_argument('--star-sj', help='Optional: STAR SJ.out.tab to convert into a splice junction track')
    parser.add_argument('--CAGE_peak', dest='cage_peak', help='Optional: CAGE peaks BED file for TSS validation track')
    parser.add_argument('--polyA_peak', dest='polya_peak', help='Optional: PolyA peaks BED file for TTS validation track')
    parser.add_argument('--refGTF', dest='ref_gtf', help='Optional: Reference GTF file for direct comparison with SQANTI3 transcripts')
    parser.add_argument('--validate-only', action='store_true', help='Validate tools and inputs only, then exit')
    parser.add_argument('--dry-run', action='store_true', help='Prepare intermediates (BED with classification) and exit before bigBed/hub generation')
    parser.add_argument('--keep-temp', action='store_true', help='Keep temporary files for debugging')
    parser.add_argument('--tables', action='store_true', help='Generate interactive HTML table reports for each category')
    parser.add_argument('--sort-by', 
                        choices=['iso_exp', 'length', 'FL', 'diff_to_TSS', 'diff_to_TTS', 
                                 'diff_to_gene_TSS', 'diff_to_gene_TTS', 'dist_to_CAGE_peak', 
                                 'dist_to_polyA_site', 'none'],
                        default='none',
                        help='Sort isoforms within each reference transcript by this metric. '
                             'Default: none (sort by genomic position only). Options: iso_exp (highest expression first), '
                             'length (longest first), FL (most full-length reads first), diff_to_TSS, diff_to_TTS, '
                             'diff_to_gene_TSS, diff_to_gene_TTS, dist_to_CAGE_peak, dist_to_polyA_site '
                             '(smallest distance first for distance metrics)')
    parser.add_argument('--no-category-tracks', action='store_true',
                        help='Only generate the main SQANTI3 track without separate tracks for each structural category')
    
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
        star_sj=args.star_sj,
        cage_peaks=args.cage_peak,
        polya_peaks=args.polya_peak,
        ref_gtf=args.ref_gtf,
        two_bit_file=args.twobit,
        validate_only=args.validate_only,
        dry_run=args.dry_run,
        sort_by=args.sort_by,
        no_category_tracks=args.no_category_tracks
    )
    converter.keep_temp = args.keep_temp
    success = converter.run()
    
    if success and args.tables:
        if generate_html_reports:
            logger.info("Generating HTML table reports...")
            reports_dir = os.path.join(args.output, "table_reports")
            try:
                generate_html_reports(args.classification, reports_dir)
            except Exception as e:
                logger.error(f"Error generating table reports: {e}")
        else:
            logger.warning("Could not import filter_isoforms.py. Skipping table generation.")
    
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
