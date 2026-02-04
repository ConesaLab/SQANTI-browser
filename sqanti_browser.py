#!/usr/bin/env python3
"""
SQANTI-browser: SQANTI3 to UCSC Genome Browser Hub Converter

Converts SQANTI3 output files (*_corrected.gtf and *_classification.txt)
to bigBed format for visualization in the UCSC Genome Browser with hub functionality.

Usage:
    python -m sqanti_browser --gtf <gtf_file> --classification <classification_file> --output <output_dir> --genome <genome>

Author: Carolina Monzo
"""

from __future__ import annotations

import argparse
import os
import sys
import subprocess
import tempfile
import shutil
import json
from pathlib import Path
from typing import Any, Optional
import pandas as pd
import numpy as np
import logging
import re
from collections import defaultdict

from src.bed_processor import BedProcessor
from src.constants import CATEGORY_ABBREV_TO_FULL
from src.validation_tracks import ValidationTrackBuilder
from src.hub_generator import HubGenerator

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Try to import the report generator
try:
    from src.filter_isoforms import generate_html_reports
except ImportError:
    generate_html_reports = None


# ============================================================================
# Custom Palette Helper Functions
# ============================================================================

def parse_color_value(color_value: list | tuple | str) -> tuple[int, int, int]:
    """Parse a color value from RGB array or hex string to RGB tuple.
    
    Args:
        color_value: Either an RGB array [R, G, B] or hex string "#RRGGBB" or "RRGGBB"
        
    Returns:
        Tuple of (R, G, B) integers
        
    Raises:
        ValueError: If color format is invalid
    """
    if isinstance(color_value, (list, tuple)):
        if len(color_value) != 3:
            raise ValueError(f"RGB array must have exactly 3 values, got {len(color_value)}")
        r, g, b = color_value
        if not all(isinstance(v, (int, float)) and 0 <= v <= 255 for v in [r, g, b]):
            raise ValueError(f"RGB values must be integers between 0-255, got {color_value}")
        return (int(r), int(g), int(b))
    elif isinstance(color_value, str):
        # Handle hex strings
        hex_str = color_value.lstrip('#')
        if len(hex_str) != 6:
            raise ValueError(f"Hex color must be 6 characters (RRGGBB), got '{color_value}'")
        try:
            r = int(hex_str[0:2], 16)
            g = int(hex_str[2:4], 16)
            b = int(hex_str[4:6], 16)
            return (r, g, b)
        except ValueError:
            raise ValueError(f"Invalid hex color string: '{color_value}'")
    else:
        raise ValueError(f"Color must be RGB array or hex string, got {type(color_value)}")


def load_custom_palette(palette_file: str | Path) -> tuple[dict[str, tuple[int, int, int]], dict[str, tuple[int, int, int]], dict[str, tuple[int, int, int]]]:
    """Load and validate a custom color palette from a JSON file.
    
    Args:
        palette_file: Path to JSON file containing palette definition
        
    Returns:
        Tuple of (standard_colors, highlight_colors) dictionaries mapping
        category names to RGB tuples
        
    Raises:
        ValueError: If JSON is invalid or contains invalid colors
        FileNotFoundError: If palette file doesn't exist
    """
    valid_categories = {
        "full-splice_match", "incomplete-splice_match", "novel_in_catalog",
        "novel_not_in_catalog", "genic", "antisense", "fusion", 
        "intergenic", "genic_intron"
    }
    
    with open(palette_file, 'r') as f:
        try:
            palette_data = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in palette file: {e}")
    
    standard_colors = {}
    highlight_colors = {}
    
    # Parse standard colors
    if "standard" in palette_data:
        for category, color_val in palette_data["standard"].items():
            if category not in valid_categories:
                logger.warning(f"Unknown category '{category}' in palette file, skipping")
                continue
            try:
                standard_colors[category] = parse_color_value(color_val)
            except ValueError as e:
                raise ValueError(f"Invalid color for standard/{category}: {e}")
    
    # Parse highlight colors (optional)
    if "highlight" in palette_data:
        for category, color_val in palette_data["highlight"].items():
            if category not in valid_categories:
                logger.warning(f"Unknown category '{category}' in highlight palette, skipping")
                continue
            try:
                highlight_colors[category] = parse_color_value(color_val)
            except ValueError as e:
                raise ValueError(f"Invalid color for highlight/{category}: {e}")
    
    # Parse validation track colors (optional)
    valid_validation = {"CAGE_peaks", "polyA_peaks", "star_sj", "reference"}
    validation_colors = {}
    if "validation_tracks" in palette_data:
        for track, color_val in palette_data["validation_tracks"].items():
            if track not in valid_validation:
                logger.warning(f"Unknown validation track '{track}' in palette, skipping")
                continue
            try:
                validation_colors[track] = parse_color_value(color_val)
            except ValueError as e:
                raise ValueError(f"Invalid color for validation_tracks/{track}: {e}")
    
    return standard_colors, highlight_colors, validation_colors


class SQANTI3ToBigBed:
    def __init__(self, gtf_file, classification_file, output_dir, genome, chrom_sizes_file=None, star_sj=None, cage_peaks=None, polya_peaks=None, ref_gtf=None, two_bit_file=None, validate_only=False, dry_run=False, sort_by='none', no_category_tracks=False, category_tracks=None, no_highlight=False, custom_palette=None):
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
        self.category_tracks = category_tracks  # None = all, or set of full category names
        self.no_highlight = no_highlight
        self.two_bit_file = two_bit_file
        self.validate_only = validate_only
        self.dry_run = dry_run
        self.keep_temp = False
        self.category_bigbeds = {}
        self.sort_by = sort_by
        
        # Initialize custom palette (will be merged with defaults)
        self.custom_standard_colors = {}
        self.custom_highlight_colors = {}
        self.custom_validation_colors = {}
        if custom_palette:
            try:
                result = load_custom_palette(custom_palette)
                self.custom_standard_colors, self.custom_highlight_colors = result[0], result[1]
                if len(result) > 2:
                    self.custom_validation_colors = result[2]
                logger.info(f"Loaded custom palette from {custom_palette}")
                if self.custom_standard_colors:
                    logger.info(f"  Custom standard colors: {list(self.custom_standard_colors.keys())}")
                if self.custom_highlight_colors:
                    logger.info(f"  Custom highlight colors: {list(self.custom_highlight_colors.keys())}")
                if self.custom_validation_colors:
                    logger.info(f"  Custom validation track colors: {list(self.custom_validation_colors.keys())}")
            except (ValueError, FileNotFoundError) as e:
                logger.error(f"Failed to load custom palette: {e}")
                sys.exit(1)
        
        # Create output directory (handles both relative and absolute paths)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Output directory: {self.output_dir.absolute()}")
        
        # Sub-components for BED processing, validation tracks, and hub generation
        self._bed_processor = BedProcessor(self)
        self._validation_builder = ValidationTrackBuilder(self)
        self._hub_generator = HubGenerator(self)
        
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
        return self._bed_processor.convert_gtf_to_genepred()

    def convert_genepred_to_bed(self, genepred_file):
        """Convert GenePred to BED format and fix malformed block arrays"""
        return self._bed_processor.convert_genepred_to_bed(genepred_file)

    def add_classification_data_to_bed(self, bed_file):
        """Add classification data to BED file using Pandas"""
        return self._bed_processor.add_classification_data_to_bed(bed_file)

    def _create_star_sj_bigbed(self, star_sj_file):
        """Convert STAR splice junctions to bigBed"""
        return self._validation_builder.create_star_sj_bigbed(star_sj_file)

    def _create_cage_bigbed(self, cage_file):
        """Convert CAGE peaks BED to bigBed"""
        return self._validation_builder.create_cage_bigbed(cage_file)

    def _create_polya_bigbed(self, polya_file):
        """Convert PolyA peaks BED to bigBed"""
        return self._validation_builder.create_polya_bigbed(polya_file)

    def _create_reference_bigbed(self, ref_gtf_file):
        """Convert reference GTF to bigBed for direct comparison"""
        return self._validation_builder.create_reference_bigbed(ref_gtf_file)

    def _generate_trix_index(self, enhanced_bed_file, genome_dir):
        """Generate Trix (.ix/.ixx) text index for fast search."""
        return self._hub_generator._generate_trix_index(enhanced_bed_file, genome_dir)

    def create_bigbed_file(self, bed_file):
        """Convert full BED to bigBed using dynamic autoSql"""
        return self._bed_processor.create_bigbed_file(bed_file)

    def create_category_bigbeds(self, main_bed_file):
        """Create separate bigBed files for each structural category"""
        return self._hub_generator.create_category_bigbeds(main_bed_file)

    def create_hub_files(self, bigbed_file):
        """Create UCSC Genome Browser hub files"""
        return self._hub_generator.create_hub_files(bigbed_file)

    def run(self) -> bool:
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
    parser.add_argument('--CAGE-peak', dest='cage_peak', help='Optional: CAGE peaks BED file for TSS validation track')
    parser.add_argument('--polyA-peak', dest='polya_peak', help='Optional: PolyA peaks BED file for TTS validation track')
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
                             'Default: none (sort by genomic position only). Options: '
                             'iso_exp (short read expression for this isoform), length (longest first), '
                             'FL (most full-length reads first), diff_to_TSS, diff_to_TTS, '
                             'diff_to_gene_TSS, diff_to_gene_TTS, dist_to_CAGE_peak, dist_to_polyA_site '
                             '(smallest distance first for distance metrics)')
    parser.add_argument('--no-category-tracks', action='store_true',
                        help='Only generate the main SQANTI3 track without separate tracks for each structural category')
    parser.add_argument('--category-tracks', type=str, metavar='LIST',
                        help='Only create tracks for these structural categories (comma-separated). '
                             'Use abbreviated names: FSM, ISM, NIC, NNC, antisense, genic_intron, genic_genomic (or genic), intergenic, fusion. '
                             'Example: --category-tracks FSM,ISM,NIC')
    parser.add_argument('--no-highlight', action='store_true',
                        help='Disable highlight coloring for top FL isoforms')
    parser.add_argument('--my-palette', type=str, metavar='JSON_FILE',
                        help='Custom color palette JSON file. Supports RGB arrays [R,G,B] or hex strings "#RRGGBB". '
                             'See ./example/example_palette.json for format. Partial specification allowed (missing categories use defaults).')
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not os.path.exists(args.gtf):
        logger.error(f"GTF file not found: {args.gtf}")
        sys.exit(1)
    
    if not os.path.exists(args.classification):
        logger.error(f"Classification file not found: {args.classification}")
        sys.exit(1)
    
    # Check if custom palette file exists
    if args.my_palette and not os.path.exists(args.my_palette):
        logger.error(f"Custom palette file not found: {args.my_palette}")
        sys.exit(1)
    
    # Parse --category-tracks (abbreviated names -> full category names)
    category_tracks = None
    if args.category_tracks:
        allowed = set()
        for abbr in (s.strip() for s in args.category_tracks.split(',') if s.strip()):
            full = CATEGORY_ABBREV_TO_FULL.get(abbr)
            if full:
                allowed.add(full)
            else:
                logger.warning(f"Unknown category abbreviation '{abbr}', skipping. Valid: FSM, ISM, NIC, NNC, antisense, genic_intron, genic_genomic, genic, intergenic, fusion")
        if allowed:
            category_tracks = allowed
            logger.info(f"Creating category tracks only for: {', '.join(sorted(allowed))}")
    
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
        no_category_tracks=args.no_category_tracks,
        category_tracks=category_tracks,
        no_highlight=args.no_highlight,
        custom_palette=args.my_palette
    )
    converter.keep_temp = args.keep_temp
    success = converter.run()
    
    if success and args.tables:
        if generate_html_reports:
            logger.info("Generating HTML table reports...")
            reports_dir = os.path.join(args.output, "table_reports")
            try:
                # Pass custom palette to HTML reports if available
                custom_palette = None
                if hasattr(converter, 'final_standard_palette'):
                    custom_palette = converter.final_standard_palette
                generate_html_reports(args.classification, reports_dir, custom_palette=custom_palette)
            except Exception as e:
                logger.error(f"Error generating table reports: {e}")
        else:
            logger.warning("Could not import src.filter_isoforms. Skipping table generation.")
    
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
