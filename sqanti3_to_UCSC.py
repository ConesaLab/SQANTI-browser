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
    def __init__(self, gtf_file, classification_file, output_dir, genome, chrom_sizes_file=None, github_repo=None, github_branch="main", enable_trix=True, star_sj=None, two_bit_file=None, validate_only=False, dry_run=False, bb12plus=False):
        self.gtf_file = gtf_file
        self.classification_file = classification_file
        self.output_dir = Path(output_dir)
        self.genome = genome
        self.chrom_sizes_file = chrom_sizes_file
        self.github_repo = github_repo
        self.github_branch = github_branch
        self.temp_dir = None
        self.enable_trix = enable_trix
        self.star_sj = star_sj
        self.star_bigbed = None
        self.two_bit_file = two_bit_file
        self.validate_only = validate_only
        self.dry_run = dry_run
        self.bb12plus = bb12plus
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
            df = pd.read_csv(self.classification_file, sep='\t')
            logger.info(f"Found {len(df)} transcripts with {len(df.columns)} columns")
            # Store column information for filter values (expects standard SQANTI3 columns)
            self.columns = df.columns.tolist()
            self.classification_data = df.set_index('isoform').to_dict('index')
            
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
    

    
    def _get_github_raw_url(self, filename):
        """Generate GitHub raw URL for a file if GitHub repo is specified"""
        if self.github_repo:
            return f"https://raw.githubusercontent.com/{self.github_repo}/{self.github_branch}/{filename}"
        else:
            return filename
    
    def copy_classification_file(self):
        """Copy the original classification file to the output directory"""
        logger.info("Copying classification file...")
        
        # Copy the original classification file
        import shutil
        classification_copy = self.output_dir / f"{self.genome}_classification.txt"
        shutil.copy2(self.classification_file, classification_copy)
        
        logger.info(f"Classification file copied: {classification_copy}")
        return classification_copy
    
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
        """Convert BED to bigBed. If bb12plus is enabled, use autoSql bed12+8 and filterByRange path"""
        logger.info("Converting BED to bigBed...")
        
        try:
            # Determine chromosome sizes source: priority --chrom-sizes > --twobit > GTF-derived
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

            # Sort BED as recommended by UCSC (LC_COLLATE=C sort -k1,1 -k2,2n)
            sorted_bed_file = os.path.join(self.temp_dir, "transcripts_with_classification.sorted.bed")
            env = os.environ.copy()
            env["LC_COLLATE"] = "C"
            sort_cmd = ['sort', '-k1,1', '-k2,2n', bed_file, '-o', sorted_bed_file]
            subprocess.run(sort_cmd, check=True, capture_output=True, text=True, env=env)
            logger.info(f"Sorted BED written to: {sorted_bed_file}")

            bigbed_file = self.output_dir / f"{self.genome}_sqanti3.bb"
            # Remember chrom sizes used for later splitting
            self.chrom_sizes_used = chrom_sizes_file
            if self.bb12plus:
                # Build a bed12+8 file using classification data columns
                bed20_path = os.path.join(self.temp_dir, 'transcripts_20col.bed')
                with open(sorted_bed_file, 'r') as inp, open(bed20_path, 'w') as outp:
                    for line in inp:
                        parts = line.rstrip('\n').split('\t')
                        if len(parts) < 12:
                            continue
                        # Compact names already hold the transcript_id
                        transcript_id = parts[3]
                        data = self.classification_data.get(transcript_id, {})
                        # pack itemRgb to uint as required by autoSql
                        try:
                            rgb_parts = parts[8].split(',')
                            if len(rgb_parts) == 3:
                                r = int(rgb_parts[0]); g = int(rgb_parts[1]); b = int(rgb_parts[2])
                                packed_rgb = (r << 16) + (g << 8) + b
                                parts[8] = str(packed_rgb)
                            else:
                                parts[8] = '0'
                        except Exception:
                            parts[8] = '0'
                        # Human-readable categorical fields
                        struct_val = data.get('structural_category','unknown')
                        subcat_val = data.get('subcategory','unknown')
                        coding_val = data.get('coding','unknown')
                        fsm = data.get('FSM_class','0')
                        # numeric fields
                        def to_int(x):
                            try:
                                return int(float(x))
                            except:
                                return 0
                        def to_float(x):
                            try:
                                return float(x)
                            except:
                                return 0.0
                        length = to_int(data.get('length',0))
                        exons = to_int(data.get('exons',0))
                        coverage = to_float(data.get('min_cov',0.0))
                        expression = to_float(data.get('iso_exp',0.0))
                        extra = [struct_val, subcat_val, coding_val, f"FSM{fsm}", str(length), str(exons), str(coverage), str(expression)]
                        outp.write('\t'.join(parts[:12] + extra) + '\n')

                # Write autoSql schema file
                as_path = os.path.join(self.temp_dir, 'sqanti3_schema.as')
                with open(as_path, 'w') as f:
                    f.write("table sqanti3Transcripts\n")
                    f.write('"SQANTI3 transcript annotations with classification data"\n')
                    f.write("(\n")
                    f.write("    string chrom;        \"Reference sequence chromosome or scaffold\"\n")
                    f.write("    uint   chromStart;   \"Start position in chromosome\"\n")
                    f.write("    uint   chromEnd;     \"End position in chromosome\"\n")
                    f.write("    string name;         \"Transcript ID\"\n")
                    f.write("    uint   score;        \"Score (0-1000)\"\n")
                    f.write("    char[1] strand;      \"+ or -\"\n")
                    f.write("    uint   thickStart;   \"Start of thick display\"\n")
                    f.write("    uint   thickEnd;     \"End of thick display\"\n")
                    f.write("    uint   itemRgb;      \"Item color packed as R*256^2+G*256+B\"\n")
                    f.write("    int    blockCount;   \"Number of blocks (exons)\"\n")
                    f.write("    int[blockCount] blockSizes;  \"Comma separated list of block sizes\"\n")
                    f.write("    int[blockCount] chromStarts; \"Start positions relative to chromStart\"\n")
                    f.write("    string structCat;    \"Structural category (e.g., full-splice_match, fusion)\"\n")
                    f.write("    string subcategory;  \"Subcategory (e.g., mono-exon, canonical)\"\n")
                    f.write("    string coding;       \"Coding status (coding, non_coding, partial_coding, pseudo)\"\n")
                    f.write("    string fsmClass;     \"FSM class (A,B,C,D)\"\n")
                    f.write("    uint   length;       \"Transcript length\"\n")
                    f.write("    uint   exons;        \"Exon count\"\n")
                    f.write("    float  coverage;     \"Minimum coverage\"\n")
                    f.write("    float  expression;   \"Isoform expression\"\n")
                    f.write(")\n")

                # Sort 20-col BED
                bed20_sorted = os.path.join(self.temp_dir, 'transcripts_20col.sorted.bed')
                env2 = os.environ.copy()
                env2["LC_COLLATE"] = "C"
                subprocess.run(['sort','-k1,1','-k2,2n',bed20_path,'-o',bed20_sorted], check=True, capture_output=True, text=True, env=env2)

                cmd = ['bedToBigBed', '-as='+as_path, '-type=bed12+8', '-extraIndex=name', bed20_sorted, chrom_sizes_file, str(bigbed_file)]
            else:
                # standard bigBed with name encoding
                cmd = ['bedToBigBed','-extraIndex=name',sorted_bed_file,chrom_sizes_file,str(bigbed_file)]
            
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"BigBed file created: {bigbed_file}")
            # Also generate per-structural-category bigBeds for hubCheck-compliant filtering via multiple tracks
            try:
                self._create_category_bigbeds(bed_file, chrom_sizes_file)
            except Exception as _e:
                logger.warning(f"Could not create per-category bigBeds: {_e}")
            return bigbed_file
            
        except Exception as e:
            logger.error(f"Error in create_bigbed_file: {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            raise e

    def _create_category_bigbeds(self, enhanced_bed_file: str, chrom_sizes_file: str) -> None:
        """Create per-structural-category bigBeds from the enhanced BED, as bed12+8 with classification fields."""
        categories = [
            'full-splice_match',
            'incomplete-splice_match',
            'novel_in_catalog',
            'novel_not_in_catalog',
            'genic',
            'antisense',
            'fusion',
            'intergenic',
            'genic_intron'
        ]
        # Create autoSql schema for bed12+8 (structCat, subcategory, coding, fsmClass, length, exons, coverage, expression)
        cat_as_path = os.path.join(self.temp_dir, 'sqanti3_cat.as')
        with open(cat_as_path, 'w') as asf:
            asf.write("table sqanti3Cat\n")
            asf.write('"SQANTI3 category tracks with classification fields"\n')
            asf.write("(\n")
            asf.write("    string chrom;        \"Reference sequence chromosome or scaffold\"\n")
            asf.write("    uint   chromStart;   \"Start position in chromosome\"\n")
            asf.write("    uint   chromEnd;     \"End position in chromosome\"\n")
            asf.write("    string name;         \"Name of item\"\n")
            asf.write("    uint   score;        \"Score (0-1000)\"\n")
            asf.write("    char[1] strand;      \"+ or -\"\n")
            asf.write("    uint   thickStart;   \"Start of thick display\"\n")
            asf.write("    uint   thickEnd;     \"End of thick display\"\n")
            asf.write("    uint   itemRgb;      \"Item color packed as R*256^2+G*256+B\"\n")
            asf.write("    int    blockCount;   \"Number of blocks\"\n")
            asf.write("    int[blockCount] blockSizes;  \"Comma separated list of block sizes\"\n")
            asf.write("    int[blockCount] chromStarts; \"Start positions relative to chromStart\"\n")
            asf.write("    string structCat;    \"Structural category\"\n")
            asf.write("    string subcategory;  \"Subcategory\"\n")
            asf.write("    string coding;       \"Coding status\"\n")
            asf.write("    string fsmClass;     \"FSM class (A,B,C,D)\"\n")
            asf.write("    uint   length;       \"Transcript length\"\n")
            asf.write("    uint   exons;        \"Exon count\"\n")
            asf.write("    float  coverage;     \"Minimum coverage (alias of min_cov)\"\n")
            asf.write("    float  expression;   \"Isoform expression (alias of iso_exp)\"\n")
            asf.write("    string associated_gene;       \"Associated gene\"\n")
            asf.write("    string associated_transcript; \"Associated transcript\"\n")
            asf.write("    uint   ref_length;            \"Reference transcript length\"\n")
            asf.write("    uint   ref_exons;             \"Reference exon count\"\n")
            asf.write("    int    diff_to_TSS;           \"Distance to TSS\"\n")
            asf.write("    int    diff_to_TTS;           \"Distance to TTS\"\n")
            asf.write("    int    diff_to_gene_TSS;      \"Distance to gene TSS\"\n")
            asf.write("    int    diff_to_gene_TTS;      \"Distance to gene TTS\"\n")
            asf.write("    string RTS_stage;             \"RTS stage\"\n")
            asf.write("    string all_canonical;         \"All junctions canonical\"\n")
            asf.write("    float  min_sample_cov;        \"Minimum sample coverage\"\n")
            asf.write("    float  min_cov;               \"Minimum coverage\"\n")
            asf.write("    uint   min_cov_pos;           \"Position of minimum coverage\"\n")
            asf.write("    float  sd_cov;                \"Coverage standard deviation\"\n")
            asf.write("    string FL;                    \"Full-length flag\"\n")
            asf.write("    uint   n_indels;              \"Number of indels\"\n")
            asf.write("    uint   n_indels_junc;         \"Number of junction indels\"\n")
            asf.write("    string bite;                  \"Bite\"\n")
            asf.write("    float  iso_exp;               \"Isoform expression\"\n")
            asf.write("    float  gene_exp;              \"Gene expression\"\n")
            asf.write("    float  ratio_exp;             \"Isoform/gene expression ratio\"\n")
            asf.write("    string FSM_class;             \"FSM class (original field)\"\n")
            asf.write("    uint   ORF_length;            \"ORF length\"\n")
            asf.write("    uint   CDS_length;            \"CDS length\"\n")
            asf.write("    uint   CDS_start;             \"CDS start\"\n")
            asf.write("    uint   CDS_end;               \"CDS end\"\n")
            asf.write("    uint   CDS_genomic_start;     \"CDS genomic start\"\n")
            asf.write("    uint   CDS_genomic_end;       \"CDS genomic end\"\n")
            asf.write("    string predicted_NMD;         \"Predicted NMD\"\n")
            asf.write("    float  perc_A_downstream_TTS; \"Percent A downstream of TTS\"\n")
            asf.write("    string seq_A_downstream_TTS;  \"Sequence downstream of TTS\"\n")
            asf.write("    int    dist_to_CAGE_peak;     \"Distance to CAGE peak\"\n")
            asf.write("    string within_CAGE_peak;      \"Within CAGE peak\"\n")
            asf.write("    int    dist_to_polyA_site;    \"Distance to polyA site\"\n")
            asf.write("    string within_polyA_site;     \"Within polyA site\"\n")
            asf.write("    string polyA_motif;           \"polyA motif\"\n")
            asf.write("    int    polyA_dist;            \"Distance to polyA motif\"\n")
            asf.write("    string polyA_motif_found;     \"polyA motif found\"\n")
            asf.write("    float  ratio_TSS;             \"TSS ratio\"\n")
            asf.write(")\n")
        # Prepare temp files per category
        cat_to_temp = {cat: os.path.join(self.temp_dir, f"cat_{cat}.bed") for cat in categories}
        # Write entries per category, appending 8 classification fields
        with open(enhanced_bed_file, 'r') as inp:
            outputs = {cat: open(path, 'w') for cat, path in cat_to_temp.items()}
            try:
                for line in inp:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 12:
                        continue
                    # recover transcript_id and struct category from encoded name if present in name
                    # name format: transcript_id|struct_cat|...
                    name = parts[3]
                    struct_cat = None
                    if '|' in name:
                        toks = name.split('|')
                        if len(toks) > 1:
                            struct_cat = toks[1]
                    # Fallback: try classification_data by transcript id
                    if not struct_cat:
                        tid = name.split('|', 1)[0]
                        data = self.classification_data.get(tid, {})
                        struct_cat = data.get('structural_category', '')
                    if struct_cat in outputs:
                        # Set category color, packed as uint for itemRgb per AS
                        try:
                            rgb = self._get_rgb_color(struct_cat) or '0,0,0'
                            r, g, b = [int(x) for x in rgb.split(',')]
                            parts[8] = str((r << 16) + (g << 8) + b)
                        except Exception:
                            parts[8] = '0'
                        # Lookup classification values
                        tid2 = parts[3].split('|', 1)[0] if '|' in parts[3] else parts[3]
                        d2 = self.classification_data.get(tid2, {})
                        def to_int(x):
                            try:
                                return int(float(x))
                            except Exception:
                                return 0
                        def to_float(x):
                            try:
                                return float(x)
                            except Exception:
                                return 0.0
                        struct_val = d2.get('structural_category', struct_cat)
                        subcat_val = d2.get('subcategory', 'unknown')
                        coding_val = d2.get('coding', 'unknown')
                        fsm_val = d2.get('FSM_class', '')
                        length_val = to_int(d2.get('length', 0))
                        try:
                            exons_val = int(float(d2.get('exons', '')))
                        except Exception:
                            try:
                                exons_val = int(parts[9])
                            except Exception:
                                exons_val = 0
                        cov_val = to_float(d2.get('min_cov', 0))
                        exp_val = to_float(d2.get('iso_exp', 0))
                        # Additional requested fields (excluding ORF_seq to keep files small)
                        extra_more = [
                            str(d2.get('associated_gene', '')),
                            str(d2.get('associated_transcript', '')),
                            str(to_int(d2.get('ref_length', 0))),
                            str(to_int(d2.get('ref_exons', 0))),
                            str(to_int(d2.get('diff_to_TSS', 0))),
                            str(to_int(d2.get('diff_to_TTS', 0))),
                            str(to_int(d2.get('diff_to_gene_TSS', 0))),
                            str(to_int(d2.get('diff_to_gene_TTS', 0))),
                            str(d2.get('RTS_stage', '')),
                            str(d2.get('all_canonical', '')),
                            str(to_float(d2.get('min_sample_cov', 0))),
                            str(cov_val),
                            str(to_int(d2.get('min_cov_pos', 0))),
                            str(to_float(d2.get('sd_cov', 0))),
                            str(d2.get('FL', '')),
                            str(to_int(d2.get('n_indels', 0))),
                            str(to_int(d2.get('n_indels_junc', 0))),
                            str(d2.get('bite', '')),
                            str(exp_val),
                            str(to_float(d2.get('gene_exp', 0))),
                            str(to_float(d2.get('ratio_exp', 0))),
                            str(d2.get('FSM_class', '')),
                            str(to_int(d2.get('ORF_length', 0))),
                            str(to_int(d2.get('CDS_length', 0))),
                            str(to_int(d2.get('CDS_start', 0))),
                            str(to_int(d2.get('CDS_end', 0))),
                            str(to_int(d2.get('CDS_genomic_start', 0))),
                            str(to_int(d2.get('CDS_genomic_end', 0))),
                            str(d2.get('predicted_NMD', '')),
                            str(to_float(d2.get('perc_A_downstream_TTS', 0))),
                            str(d2.get('seq_A_downstream_TTS', '')),
                            str(to_int(d2.get('dist_to_CAGE_peak', 0))),
                            str(d2.get('within_CAGE_peak', '')),
                            str(to_int(d2.get('dist_to_polyA_site', 0))),
                            str(d2.get('within_polyA_site', '')),
                            str(d2.get('polyA_motif', '')),
                            str(to_int(d2.get('polyA_dist', 0))),
                            str(d2.get('polyA_motif_found', '')),
                            str(to_float(d2.get('ratio_TSS', 0))),
                        ]
                        extra = [
                            str(struct_val), str(subcat_val), str(coding_val), str(fsm_val),
                            str(length_val), str(exons_val), str(cov_val), str(exp_val)
                        ] + extra_more
                        outputs[struct_cat].write('\t'.join(parts[:12] + extra) + '\n')
            finally:
                for fh in outputs.values():
                    fh.close()

        # Convert each category BED to sorted bigBed (bed12+8)
        for cat, bed_path in cat_to_temp.items():
            # Skip empty files
            if not os.path.exists(bed_path) or os.path.getsize(bed_path) == 0:
                continue
            sorted_path = os.path.join(self.temp_dir, f"cat_{cat}.sorted.bed")
            env = os.environ.copy()
            env["LC_COLLATE"] = "C"
            subprocess.run(['sort','-k1,1','-k2,2n', bed_path, '-o', sorted_path], check=True, capture_output=True, text=True, env=env)
            out_bb = self.output_dir / f"{self.genome}_sqanti3_{cat}.bb"
            # bed12 + 47 extra fields (see sqanti3_cat.as)
            subprocess.run(['bedToBigBed', '-as='+cat_as_path, '-type=bed12+47', sorted_path, chrom_sizes_file, str(out_bb)], check=True, capture_output=True, text=True)
            self.category_bigbeds[cat] = out_bb

    def add_classification_data_to_bed(self, bed_file):
        """Add classification data to name field for smart filtering in standard bigBed"""
        logger.info("Adding classification data to BED file...")
        
        enhanced_bed_file = os.path.join(self.temp_dir, "transcripts_with_classification.bed")
        processed_count = 0
        missing_count = 0
        
        # Define structural category encoding (numeric codes)
        structural_categories = [
            'full-splice_match',
            'incomplete-splice_match', 
            'novel_in_catalog',
            'novel_not_in_catalog',
            'genic',
            'antisense',
            'fusion',
            'intergenic',
            'genic_intron'
        ]
        
        # Define subcategory encoding (numeric codes)
        subcategories = [
            'mono-exon',
            'multi-exon',
            'novel',
            'known',
            'canonical',
            'non-canonical'
        ]
        
        # Define coding status encoding (numeric codes)
        coding_statuses = [
            'coding',
            'non_coding',
            'partial_coding',
            'pseudo'
        ]
        
        with open(bed_file, 'r') as infile, open(enhanced_bed_file, 'w') as outfile:
            for line in infile:
                parts = line.strip().split('\t')
                if len(parts) >= 12:  # BED12 format
                    transcript_id = parts[3]
                    
                    # Look up classification data
                    if transcript_id in self.classification_data:
                        data = self.classification_data[transcript_id]
                        structural_category = data.get('structural_category', '')
                        
                        # Get RGB color for structural category and write as R,G,B string (UCSC expected format)
                        rgb_color = self._get_rgb_color(structural_category)
                        parts[8] = rgb_color if rgb_color else "0"
                        
                        # Encode classification data in name field for comprehensive filtering
                        # Format: transcript_id|struct_cat|subcat|coding|fsm_class|length|exons|coverage|expression
                        struct_code = 0  # Default
                        for i, cat in enumerate(structural_categories):
                            if structural_category.lower() == cat.lower():
                                struct_code = i
                                break
                        
                        subcategory_val = data.get('subcategory', '')
                        subcat_code = 0  # Default
                        for i, subcat in enumerate(subcategories):
                            if subcategory_val.lower() == subcat.lower():
                                subcat_code = i
                                break
                        
                        coding_val = data.get('coding', '')
                        coding_code = 0  # Default
                        for i, code in enumerate(coding_statuses):
                            if coding_val.lower() == code.lower():
                                coding_code = i
                                break
                        
                        fsm_class_val = data.get('FSM_class', '')
                        fsm_code = 0  # Default
                        if fsm_class_val in ['A', 'B', 'C', 'D']:
                            fsm_code = ord(fsm_class_val) - ord('A') + 1  # 1=A, 2=B, 3=C, 4=D
                        
                        length = data.get('length', '0')
                        try:
                            length_val = int(float(length)) if length and length != 'NA' else 0
                        except:
                            length_val = 0
                        
                        exons = data.get('exons', '0')
                        try:
                            exons_val = int(float(exons)) if exons and exons != 'NA' else 0
                        except:
                            exons_val = 0
                        
                        min_cov = data.get('min_cov', '0')
                        try:
                            cov_val = float(min_cov) if min_cov and min_cov != 'NA' else 0.0
                        except:
                            cov_val = 0.0
                        
                        iso_exp = data.get('iso_exp', '0')
                        try:
                            exp_val = float(iso_exp) if iso_exp and iso_exp != 'NA' else 0.0
                        except:
                            exp_val = 0.0
                        
                        # Set score field to exon count for score-based filtering in UCSC UI
                        try:
                            parts[4] = str(min(max(exons_val, 0), 1000))
                        except Exception:
                            parts[4] = '0'
                        
                        # Encode ALL classification data in name field for comprehensive filtering
                        # Format: transcript_id|struct_cat|subcat|coding|fsm|length|exons|coverage|expression
                        struct_cat_name = structural_categories[struct_code] if struct_code < len(structural_categories) else "unknown"
                        subcat_name = subcategories[subcat_code] if subcat_code < len(subcategories) else "unknown"
                        coding_name = coding_statuses[coding_code] if coding_code < len(coding_statuses) else "unknown"
                        
                        # Use compact name (transcript ID only); rich text goes to Trix
                        parts[3] = transcript_id
                        

                        
                        processed_count += 1
                    else:
                        # Default values for transcripts not in classification
                        parts[8] = "16777215"  # White
                        parts[4] = '0'
                        
                        # Use compact name even when classification is missing
                        parts[3] = transcript_id
                        
                        missing_count += 1
                    
                    outfile.write('\t'.join(parts) + '\n')
        
        logger.info(f"Classification data added: {processed_count} processed, {missing_count} missing from classification")
        return enhanced_bed_file

    def _generate_trix_index(self, enhanced_bed_file, genome_dir):
        """Generate Trix (.ix/.ixx) text index with rich descriptions for fast search"""
        ixixx_path = shutil.which('ixIxx')
        if not ixixx_path:
            logger.warning("ixIxx not found; skipping Trix index generation")
            return False

        try:
            trix_input = os.path.join(self.temp_dir, 'trix_input.txt')
            with open(enhanced_bed_file, 'r') as bed_in, open(trix_input, 'w') as t_out:
                for line in bed_in:
                    parts = line.rstrip('\t\n').split('\t')
                    if len(parts) >= 12:
                        tid = parts[3]
                        d = self.classification_data.get(tid, {})
                        fields = {
                            'structural_category': d.get('structural_category', 'unknown'),
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
                            'ratio_TSS': d.get('ratio_TSS', ''),
                        }
                        desc = f"{tid} " + ' '.join([f"{k} {v}" for k, v in fields.items() if v not in (None, '', 'NA')])
                        synonyms = ' '.join([str(v) for v in fields.values() if v not in (None, '', 'NA')])
                        t_out.write(f"{tid}\t{desc}\t{synonyms}\n")

            ix_path = os.path.join(genome_dir, 'trix.ix')
            ixx_path = os.path.join(genome_dir, 'trix.ixx')
            cmd = ['ixIxx', trix_input, ix_path, ixx_path]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"Trix index generated: {ix_path}, {ixx_path}")
            return True
        except Exception as e:
            logger.warning(f"Failed to generate Trix index: {e}")
            return False

    def _create_star_sj_bigbed(self, sj_file):
        """Convert STAR SJ.out.tab to bigBed splice junction track"""
        try:
            if not os.path.exists(sj_file):
                logger.error(f"STAR SJ file not found: {sj_file}")
                return None

            # Determine chromosome sizes source
            if self.chrom_sizes_file and os.path.exists(self.chrom_sizes_file):
                chrom_sizes_file = self.chrom_sizes_file
            else:
                chrom_sizes_file = self.extract_chrom_sizes()
                if not chrom_sizes_file:
                    raise Exception("Could not determine chromosome sizes for STAR junctions")

            # Create BED6 from SJ.out.tab
            bed6_path = os.path.join(self.temp_dir, 'SJ.out.bed')
            strand_map = {'0':'.', '1':'+', '2':'-'}
            with open(sj_file, 'r') as sj_in, open(bed6_path, 'w') as bed_out:
                for idx, line in enumerate(sj_in, 1):
                    parts = line.strip().split('\t')
                    if len(parts) < 4:
                        continue
                    chrom = parts[0]
                    start = int(parts[1]) - 1
                    end = int(parts[2])
                    strand_code = parts[3]
                    strand = strand_map.get(strand_code, '.')
                    name = f"j_{idx}"
                    score = '1000'
                    bed_out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

            # Sort BED
            sorted_bed_path = os.path.join(self.temp_dir, 'SJ.out.sorted.bed')
            env = os.environ.copy()
            env["LC_COLLATE"] = "C"
            sort_cmd = ['sort', '-k1,1', '-k2,2n', bed6_path, '-o', sorted_bed_path]
            subprocess.run(sort_cmd, check=True, capture_output=True, text=True, env=env)

            # Create bigBed (type bed6)
            star_bb = self.output_dir / f"{self.genome}_star_sj.bb"
            cmd = ['bedToBigBed', sorted_bed_path, chrom_sizes_file, str(star_bb)]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"STAR junctions bigBed created: {star_bb}")
            return star_bb
        except Exception as e:
            logger.error(f"Failed to create STAR junctions bigBed: {e}")
            return None

    def create_autosql_schema(self):
        """Create autoSql schema file for bigBed 12 + format"""
        logger.info("Creating autoSql schema for bigBed 12 + format...")
        
        try:
            as_file = os.path.join(self.temp_dir, "sqanti3_schema.as")
            
            with open(as_file, 'w') as f:
                f.write("table sqanti3Transcripts\n")
                f.write('"SQANTI3 transcript annotations with classification data"\n')
                f.write("(\n")
                # Standard BED12 fields
                f.write("    string chrom;       \"Reference sequence chromosome or scaffold\"\n")
                f.write("    uint   chromStart;  \"Start position in chromosome\"\n")
                f.write("    uint   chromEnd;    \"End position in chromosome\"\n")
                f.write("    string name;        \"Name of transcript\"\n")
                f.write("    uint   score;       \"Score (0-1000)\"\n")
                f.write("    char[1] strand;     \"+ or - for strand\"\n")
                f.write("    uint   thickStart;  \"Start of where display should be thick (start codon)\"\n")
                f.write("    uint   thickEnd;    \"End of where display should be thick (stop codon)\"\n")
                f.write("    uint   itemRgb;     \"RGB color value (packed)\"\n")
                f.write("    int    blockCount;  \"Number of blocks\"\n")
                f.write("    int[blockCount] blockSizes; \"Comma separated list of block sizes\"\n")
                f.write("    int[blockCount] chromStarts; \"Start positions relative to chromStart\"\n")
                
                # Custom fields for classification data (using UCSC-required field names)
                f.write("    uint   expCount;    \"Structural category code (0-8)\"\n")
                f.write("    uint   isSizeLink;  \"Subcategory code (0-5)\"\n")
                f.write("    uint   expIds;      \"Coding status code (0-3)\"\n")
                f.write("    uint   expScores;   \"FSM class code (0-4)\"\n")
                f.write("    uint   length;      \"Transcript length\"\n")
                f.write("    uint   exons;       \"Number of exons\"\n")
                f.write("    float  coverage;    \"Minimum coverage\"\n")
                f.write("    float  expression;  \"Isoform expression\"\n")
                f.write(")\n")
            
            logger.info(f"AutoSql schema created: {as_file}")
            return as_file
            
        except Exception as e:
            logger.error(f"Error creating autoSql schema: {e}")
            return None

    def create_hub_files(self, bigbed_file):
        """Create UCSC Genome Browser hub files"""
        logger.info("Creating hub files...")
        
        # Extract the output directory name for hub naming
        # Always use just the final directory name for clean, simple hub names
        output_dir_name = self.output_dir.name
        hub_name = f"{output_dir_name}_{self.genome}_SQANTI3_Hub"
        
        # Copy the original classification file
        classification_copy = self.copy_classification_file()
        
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
            f.write(f"descriptionUrl {self._get_github_raw_url(f'{hub_name}.html')}\n")
        
        # Create genomes.txt
        genomes_file = self.output_dir / "genomes.txt"
        with open(genomes_file, 'w', newline='\n') as f:
            f.write(f"genome {self.genome}\n")
            f.write(f"trackDb {self.genome}/trackDb.txt\n")
            # reference groups file for multi-track organization
            f.write(f"groups {self.genome}/groups.txt\n")
        
        # Create genome-specific directory and trackDb.txt
        genome_dir = self.output_dir / self.genome
        genome_dir.mkdir(exist_ok=True)
        
        # Create groups.txt for organizing tracks
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
<html lang=\"en\">
<head>
    <meta charset=\"UTF-8\">
    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">
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
        {('<li>Use per-field filters and ranges (length, exons, coverage, expression) with 12+ schema.</li>' if self.bb12plus else '')}
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
<html lang=\"en\">
<head>
    <meta charset=\"UTF-8\">
    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">
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
        
        trackdb_file = genome_dir / "trackDb.txt"
        with open(trackdb_file, 'w', newline='\n') as f:
            f.write(f"track {self.genome}_sqanti3\n")
            f.write(f"bigDataUrl {self._get_github_raw_url(f'{self.genome}_sqanti3.bb')}\n")
            f.write(f"shortLabel SQANTI3 Transcripts\n")
            f.write(f"longLabel SQANTI3 Transcriptome Analysis Results\n")
            if self.bb12plus:
                f.write(f"type bigBed 12 + 8\n")
            else:
                f.write(f"type bigBed\n")
            f.write(f"visibility dense\n")
            f.write(f"group transcripts\n")
            f.write(f"itemRgb on\n")
            f.write(f"color 107,174,214\n")
            f.write(f"priority 1\n")
            # Note: hubCheck warns on useScore; omit score filter directives
            f.write(f"html {self._get_github_raw_url(transcripts_html_name)}\n")
            # Enable search on name; fielded filters are added below when 12+ is enabled
            f.write(f"searchIndex name\n")
            # Add exon range control bound to built-in BED12 field blockCount
            f.write(f"filterByRange.blockCount 0:100\n")
            f.write(f"filterLabel.blockCount Number of exons\n")
            
            # Add Trix search index if present
            trix_ix = genome_dir / 'trix.ix'
            if os.path.exists(trix_ix):
                f.write(f"searchTrix trix.ix\n")

            # Append STAR junctions track if generated
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

            # Append per-structural-category tracks for hubCheck-compliant filtering via toggles
            if self.category_bigbeds:
                for cat, bb_path in self.category_bigbeds.items():
                    # Create per-category HTML description
                    cat_html_name = f"{self.genome}_sqanti3_{cat}.html"
                    try:
                        with open(self.output_dir / cat_html_name, 'w') as ch:
                            ch.write(f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
    <meta charset=\"UTF-8\">
    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">
    <title>{self.genome} SQANTI3 {cat}</title>
    <style>
        body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }}
        h1 {{ color: #333; }}
        ul {{ list-style-type: disc; margin-left: 20px; }}
    </style>
    </head>
<body>
    <h1>SQANTI3 {cat} ({self.genome})</h1>
    <p>This track shows only transcripts in the <strong>{cat}</strong> structural category.</p>
    <p>Colors follow the SQANTI3 palette. Use the track display settings to adjust visibility and score (exon count) range.</p>
</body>
</html>""")
                    except Exception:
                        pass
                    f.write("\n")
                    track_name = f"{self.genome}_sqanti3_{cat}"
                    short = f"SQANTI3 {cat}"
                    f.write(f"track {track_name}\n")
                    f.write(f"bigDataUrl {self._get_github_raw_url(bb_path.name)}\n")
                    f.write(f"shortLabel {short}\n")
                    f.write(f"longLabel SQANTI3 {cat} transcripts\n")
                    # Declare the full number of extra fields added to category tracks
                    f.write(f"type bigBed 12 + 39\n")
                    f.write(f"visibility hide\n")
                    f.write(f"group transcripts\n")
                    f.write(f"itemRgb on\n")
                    # Omit useScore/scoreFilter directives for hubCheck compatibility
                    f.write(f"priority 3\n")
                    f.write(f"html {self._get_github_raw_url(cat_html_name)}\n")
                    # Exon range label via blockCount
                    f.write(f"filterByRange.blockCount 0:100\n")
                    f.write(f"filterLabel.blockCount Number of exons\n")
        
        # Create simple HTML description
        html_file = self.output_dir / f"{hub_name}.html"
        with open(html_file, 'w') as f:
            f.write(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{hub_name} - SQANTI3 Transcriptome Analysis Results</title>
    <style>
        body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }}
        h1 {{ color: #333; }}
        h2 {{ color: #555; }}
        ul {{ list-style-type: disc; margin-left: 20px; }}
        code {{ background-color: #f4f4f4; padding: 2px 4px; border-radius: 4px; }}
        .info-box {{ background-color: #e7f3ff; padding: 15px; border-radius: 6px; margin: 20px 0; }}
    </style>
</head>
<body>
    <h1>{hub_name}</h1>
    <p>This hub displays SQANTI3 transcriptome analysis results for the {self.genome} genome assembly.</p>
    
    <div class="info-box">
        <h3>🔍 Advanced Filtering with Smart Name Encoding</h3>
        <p><strong>This hub uses standard bigBed format with comprehensive data encoding in the name field:</strong></p>
        <ul>
            <li><strong>Structural Category:</strong> Filter by FSM, ISM, NIC, NNC, genic, antisense, fusion, intergenic, genic_intron</li>
            <li><strong>Subcategory:</strong> Filter by mono-exon, multi-exon, novel, known, canonical, non-canonical</li>
            <li><strong>Coding Status:</strong> Filter by coding, non_coding, partial_coding, pseudo</li>
            <li><strong>FSM Class:</strong> Filter by A, B, C, D</li>
            <li><strong>Length:</strong> Range-based filtering for transcript length (L value)</li>
            <li><strong>Coverage:</strong> Range-based filtering for minimum coverage (C value)</li>
            <li><strong>Expression:</strong> Range-based filtering for isoform expression (X value)</li>
        </ul>
        <p><em>Use the search box to filter by classification data (e.g., "intergenic", "L>500", "coding")</em></p>
    </div>
    
    <div class="info-box">
        <h3>🔍 How the Smart Filtering Works</h3>
        <p><strong>All classification data is encoded in the transcript name field:</strong></p>
        <p><code>PB.118859.1|intergenic|mono-exon|non_coding|FSM1|L320|E1|Cnan|X80.0265</code></p>
        <ul>
            <li><strong>Format:</strong> transcript_id|structural_category|subcategory|coding_status|FSM_class|Llength|Eexons|Ccoverage|Xexpression</li>
            <li><strong>Search Examples:</strong></li>
            <ul>
                <li><code>intergenic</code> - Find all intergenic transcripts</li>
                <li><code>L>500</code> - Find transcripts longer than 500bp</li>
                <li><code>E>2</code> - Find multi-exon transcripts</li>
                <li><code>coding</code> - Find all coding transcripts</li>
                <li><code>FSM1</code> - Find FSM class A transcripts</li>
            </ul>
        </ul>
        <p><em>This approach provides comprehensive filtering without the complexity of bigBed 12+8 format!</em></p>
    </div>
    
    <h2>🎨 Color Legend</h2>
    <ul>
        <li><span style="color: #6BAED6;"><strong>Full-splice Match (FSM):</strong></span> Transcripts that perfectly match known genes</li>
        <li><span style="color: #FC8D59;"><strong>Incomplete-splice Match (ISM):</strong></span> Transcripts that partially match known genes</li>
        <li><span style="color: #78C679;"><strong>Novel In Catalog (NIC):</strong></span> Novel transcripts in known gene loci</li>
        <li><span style="color: #EE6A50;"><strong>Novel Not In Catalog (NNC):</strong></span> Novel transcripts in novel loci</li>
        <li><span style="color: #969696;"><strong>Genic:</strong></span> Transcripts in gene regions</li>
        <li><span style="color: #66C2A4;"><strong>Antisense:</strong></span> Antisense transcripts</li>
        <li><span style="color: #6BAED6;"><strong>Fusion:</strong></span> Fusion transcripts</li>
        <li><span style="color: #E9967A;"><strong>Intergenic:</strong></span> Intergenic transcripts</li>
        <li><span style="color: #41B6C4;"><strong>Genic Intron:</strong></span> Intronic transcripts</li>
    </ul>
    
    <h2>📊 Classification Data</h2>
    <p>The original SQANTI3 classification file ({classification_copy.name}) is included for reference and detailed analysis.</p>
    <p>This file contains all the detailed classification information for each transcript.</p>
    
    <h2>🚀 Usage Instructions</h2>
    <p>To use this hub in the UCSC Genome Browser:</p>
    <ol>
        <li>Upload all files in this directory to a web-accessible location</li>
        <li>In the UCSC Genome Browser, go to <strong>My Data → Track Hubs</strong></li>
        <li>Enter the URL to your <code>hub.txt</code> file</li>
        <li>Select the appropriate genome assembly</li>
        <li>The SQANTI3 tracks will appear in your track list</li>
        <li>Right-click on the track and select "Filter" to access individual field filters</li>
    </ol>
    
    <h2>📧 Contact Information</h2>
    <p><strong>Note:</strong> The hub.txt file includes a GitHub-specific email address. If you need to use a different email:</p>
    <ul>
        <li>Edit the <code>hub.txt</code> file</li>
        <li>Change the <code>email</code> line to your preferred email address</li>
        <li>Re-upload the updated file</li>
    </ul>
</body>
</html>""")
        
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

            # Optionally generate Trix index (after name encoding)
            if self.enable_trix:
                genome_dir = self.output_dir / self.genome
                genome_dir.mkdir(exist_ok=True)
                self._generate_trix_index(bed_file, genome_dir)

            # Create bigBed file
            bigbed_file = self.create_bigbed_file(bed_file)
            if not bigbed_file:
                return False
            
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
            # Only cleanup on success
            pass

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
    # Trix enabled by default; allow disabling with --no-trix
    parser.add_argument('--enable-trix', dest='enable_trix', action='store_true', default=True, help='Generate Trix (.ix/.ixx) text index for fast search (default: on)')
    parser.add_argument('--no-trix', dest='enable_trix', action='store_false', help='Disable Trix index generation')
    parser.add_argument('--star-sj', help='Optional: STAR SJ.out.tab to convert into a splice junction track')
    parser.add_argument('--bb12plus', action='store_true', help='Use bigBed 12+ schema with autoSql and filterByRange')
    parser.add_argument('--validate-only', action='store_true', help='Validate tools and inputs only, then exit')
    parser.add_argument('--dry-run', action='store_true', help='Prepare intermediates (BED with classification) and exit before bigBed/hub generation')
    
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
        enable_trix=args.enable_trix,
        star_sj=args.star_sj,
        two_bit_file=args.twobit,
        validate_only=args.validate_only,
        dry_run=args.dry_run,
        bb12plus=args.bb12plus
    )
    success = converter.run()
    
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
