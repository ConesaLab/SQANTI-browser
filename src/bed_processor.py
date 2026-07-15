"""
BED processing: GTF → GenePred → BED, classification merge, sorting, bigBed creation.
"""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path

from typing import Any, Optional

import numpy as np
import pandas as pd

from src.constants import (
    BED12_BASE_COLS,
    DEFAULT_STANDARD_PALETTE,
    DEFAULT_HIGHLIGHT_PALETTE,
)
from src.utils import darken_color

logger = logging.getLogger(__name__)


class BedProcessor:
    """Handles GTF→BED conversion, classification merge, sorting, and main bigBed creation."""

    def __init__(self, converter: Any) -> None:
        self.converter = converter

    def convert_gtf_to_genepred(self) -> Optional[str]:
        """Convert GTF to GenePred format"""
        logger.info("Converting GTF to GenePred format...")
        genepred_file = os.path.join(self.converter.temp_dir, "transcripts.genepred")
        cmd = [
            'gtfToGenePred',
            '-genePredExt',
            '-allErrors',
            '-ignoreGroupsWithoutExons',
            self.converter.gtf_file,
            genepred_file
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info("GTF to GenePred conversion completed")
            return genepred_file
        except subprocess.CalledProcessError as e:
            logger.error(f"GTF to GenePred conversion failed: {e}")
            if e.stderr:
                logger.error(f"Error details: {e.stderr}")
            return None

    def convert_genepred_to_bed(self, genepred_file: str) -> Optional[str]:
        """Convert GenePred to BED format and fix malformed block arrays"""
        logger.info("Converting GenePred to BED format...")
        try:
            temp_bed_file = os.path.join(self.converter.temp_dir, "transcripts_temp.bed")
            cmd = ['genePredToBed', genepred_file, temp_bed_file]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            bed_file = os.path.join(self.converter.temp_dir, "transcripts.bed")
            with open(temp_bed_file, 'r') as infile, open(bed_file, 'w') as outfile:
                for line in infile:
                    parts = line.strip().split('\t')
                    if len(parts) >= 12:
                        if parts[10].endswith(','):
                            parts[10] = parts[10].rstrip(',')
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
            f.write('    uint   reserved;     "Item color (itemRgb), packed as R*256^2+G*256+B"\n')
            f.write('    int    blockCount;   "Number of blocks (exons)"\n')
            f.write('    int[blockCount] blockSizes;  "Comma separated list of block sizes"\n')
            f.write('    int[blockCount] chromStarts; "Start positions relative to chromStart"\n')
            for col in extra_cols:
                safe_col = col.replace('.', '_').replace(' ', '_').replace('-', '_').replace('/', '_')
                f.write(f'    lstring {safe_col}; "{col}"\n')
            f.write(")\n")

    def add_classification_data_to_bed(self, bed_file: str) -> Optional[str]:
        """Add classification data to BED file using Pandas"""
        logger.info("Adding classification data to BED file (Vectorized)...")
        output_file = os.path.join(self.converter.temp_dir, "transcripts_full.bed")
        try:
            bed_df = pd.read_csv(
                bed_file,
                sep='\t',
                names=list(BED12_BASE_COLS),
                dtype={'chrom': 'string', 'name': 'string', 'strand': 'category'}
            )
            if not hasattr(self.converter, 'classification_df'):
                logger.error("Classification data not loaded")
                return None

            exclude_cols = set(BED12_BASE_COLS)
            exclude_cols.add('ORF_seq')
            self.converter.extra_cols = [
                c for c in self.converter.classification_df.columns
                if c not in exclude_cols and c != 'name'
            ]
            class_subset = self.converter.classification_df[['name'] + self.converter.extra_cols].copy()
            merged = bed_df.merge(class_subset, on='name', how='left')

            categorical_defaults = {
                'structural_category': 'NA',
                'subcategory': 'NA',
                'coding': 'NA',
                'FSM_class': 'NA'
            }
            for col, default_val in categorical_defaults.items():
                if col in merged.columns:
                    if isinstance(merged[col].dtype, pd.CategoricalDtype):
                        if default_val not in merged[col].cat.categories:
                            merged[col] = merged[col].cat.add_categories([default_val])
                    merged[col] = merged[col].fillna(default_val)
            for col in self.converter.extra_cols:
                if col in merged.columns:
                    merged[col] = merged[col].fillna('NA')

            def pack_rgb(r, g, b):
                return r * 65536 + g * 256 + b

            final_standard = {**DEFAULT_STANDARD_PALETTE}
            if self.converter.custom_standard_colors:
                final_standard.update(self.converter.custom_standard_colors)
            final_highlight = {**DEFAULT_HIGHLIGHT_PALETTE}
            if self.converter.custom_highlight_colors:
                final_highlight.update(self.converter.custom_highlight_colors)
            elif self.converter.custom_standard_colors:
                for cat, rgb in self.converter.custom_standard_colors.items():
                    if cat not in self.converter.custom_highlight_colors:
                        final_highlight[cat] = darken_color(rgb)

            self.converter.final_standard_palette = final_standard
            self.converter.final_highlight_palette = final_highlight
            cat_palette = {cat: pack_rgb(*rgb) for cat, rgb in final_standard.items()}
            highlight_palette = {cat: pack_rgb(*rgb) for cat, rgb in final_highlight.items()}
            default_color = pack_rgb(200, 200, 200)

            if 'structural_category' in merged.columns:
                merged['itemRgb'] = merged['structural_category'].map(cat_palette).fillna(default_color).astype(int)
            else:
                merged['itemRgb'] = default_color

            if (not self.converter.no_highlight) and all(col in merged.columns for col in ['associated_transcript', 'structural_category', 'FL']):
                fl_numeric = pd.to_numeric(merged['FL'], errors='coerce')
                gene_based_categories = {
                    "novel_in_catalog", "novel_not_in_catalog", "genic", "antisense",
                    "fusion", "intergenic", "genic_intron",
                }
                group_key = merged['associated_transcript']
                if 'associated_gene' in merged.columns:
                    use_gene = merged['structural_category'].isin(gene_based_categories)
                    group_key = group_key.where(~use_gene, merged['associated_gene'])
                valid_groups = group_key.notna() & merged['structural_category'].notna()
                if valid_groups.any():
                    max_fl = fl_numeric.groupby([group_key, merged['structural_category']]).transform('max')
                    top_mask = valid_groups & fl_numeric.notna() & max_fl.notna() & (fl_numeric == max_fl)
                    if top_mask.any():
                        merged.loc[top_mask, 'itemRgb'] = (
                            merged.loc[top_mask, 'structural_category']
                            .map(highlight_palette)
                            .fillna(merged.loc[top_mask, 'itemRgb'])
                            .astype(int)
                        )

            if 'FSM_class' in merged.columns:
                merged['FSM_class'] = merged['FSM_class'].astype(str)
                mask = merged['FSM_class'] != 'NA'
                merged.loc[mask, 'FSM_class'] = 'FSM' + merged.loc[mask, 'FSM_class']

            final_cols = list(BED12_BASE_COLS) + self.converter.extra_cols
            merged[final_cols].to_csv(output_file, sep='\t', header=False, index=False, quoting=3)
            logger.info(f"Created enhanced BED file with {len(self.converter.extra_cols)} extra columns: {output_file}")
            return output_file
        except Exception as e:
            logger.error(f"Error in vectorized BED processing: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None

    def sort_bed_for_visualization(self, df: pd.DataFrame, sort_by: Optional[str] = None) -> pd.DataFrame:
        """Sort BED DataFrame for optimal visualization in UCSC Genome Browser."""
        df = df.copy()
        df['_original_order'] = np.arange(len(df))
        if sort_by is None:
            sort_by = self.converter.sort_by
        # basic = genomic position + pipeline order for ties; none = genomic position + GTF file order for ties
        use_gtf_order = sort_by == 'none'
        if sort_by in ('basic', 'none'):
            sort_by = None
        has_associated_transcript = 'associated_transcript' in df.columns
        if not has_associated_transcript:
            logger.warning("Column 'associated_transcript' not found. Isoforms will not be grouped by reference transcript.")
        original_sort_by = sort_by
        if sort_by is not None:
            if sort_by not in df.columns:
                logger.warning(f"Sort column '{sort_by}' not found in data.")
                sort_by = None
            elif df[sort_by].isna().all() or (df[sort_by].astype(str) == 'NA').all():
                logger.warning(f"Column '{sort_by}' contains only NA values.")
                sort_by = None
            else:
                na_count = df[sort_by].isna().sum() + (df[sort_by].astype(str) == 'NA').sum()
                if na_count > 0:
                    logger.warning(f"Column '{sort_by}' has {na_count} NA values out of {len(df)} rows. NA values will be sorted last.")
        if sort_by is None and original_sort_by is not None and original_sort_by not in ('basic', 'none'):
            logger.warning(f"Could not sort by '{original_sort_by}'. Falling back to genomic position order.")

        if sort_by:
            sort_col_numeric = f'_sort_key_{sort_by}'
            df[sort_col_numeric] = pd.to_numeric(df[sort_by], errors='coerce')
        sort_cols = ['chrom', 'chromStart']
        ascending = [True, True]
        use_transcript_grouping = has_associated_transcript and sort_by is not None
        if use_transcript_grouping:
            sort_cols.append('associated_transcript')
            ascending.append(True)
        if sort_by:
            sort_cols.append(sort_col_numeric)
            if sort_by in ['diff_to_TSS', 'diff_to_TTS', 'diff_to_gene_TSS', 'diff_to_gene_TTS',
                           'dist_to_CAGE_peak', 'dist_to_polyA_site']:
                ascending.append(True)
            else:
                ascending.append(False)
        else:
            gtf_order = getattr(self.converter, 'gtf_transcript_order', None) if use_gtf_order else None
            if use_gtf_order and gtf_order:
                max_ord = len(gtf_order)
                df['_gtf_order'] = df['name'].map(gtf_order).fillna(max_ord).astype(int)
                sort_cols.append('_gtf_order')
                ascending.append(True)
            else:
                if use_gtf_order and not gtf_order:
                    logger.warning("GTF transcript order not available; using pipeline order for ties")
                sort_cols.append('_original_order')
                ascending.append(True)

        df_sorted = df.sort_values(
            by=sort_cols,
            ascending=ascending,
            na_position='last',
            key=lambda col: col if col.name != 'chrom' else col.astype(str)
        )
        if sort_by:
            df_sorted = df_sorted.drop(columns=[sort_col_numeric])
        df_sorted = df_sorted.drop(columns=['_original_order'])
        if '_gtf_order' in df_sorted.columns:
            df_sorted = df_sorted.drop(columns=['_gtf_order'])
        if sort_by:
            sort_info = f"by {sort_by}"
        elif use_gtf_order:
            sort_info = "by genomic position, GTF file order for ties"
        else:
            sort_info = "by genomic position, pipeline order for ties (basic)"
        group_info = "grouped by reference transcript, " if use_transcript_grouping else ""
        logger.info(f"Sorted {len(df_sorted)} transcripts: {group_info}{sort_info}")
        return df_sorted

    def create_bigbed_file(self, bed_file: str) -> Optional[Path]:
        """Convert full BED to bigBed using dynamic autoSql"""
        logger.info("Converting BED to bigBed...")
        result = None
        try:
            if self.converter.chrom_sizes_file and os.path.exists(self.converter.chrom_sizes_file):
                chrom_sizes_file = self.converter.chrom_sizes_file
                logger.info(f"Using provided chrom.sizes file: {chrom_sizes_file}")
            elif self.converter.two_bit_file:
                chrom_sizes_file = self.converter.extract_chrom_sizes_from_twobit()
                if not chrom_sizes_file:
                    raise Exception("twoBitInfo failed to generate chrom.sizes")
            else:
                chrom_sizes_file = self.converter.extract_chrom_sizes()
                if not chrom_sizes_file:
                    raise Exception("Could not determine chromosome sizes")

            bed_cols = list(BED12_BASE_COLS) + self.converter.extra_cols
            df = pd.read_csv(bed_file, sep='\t', names=bed_cols, dtype={'chrom': 'string'})
            df_sorted = self.sort_bed_for_visualization(df, self.converter.sort_by)
            sorted_bed_file = os.path.join(self.converter.temp_dir, "transcripts_full.sorted.bed")
            df_sorted.to_csv(sorted_bed_file, sep='\t', header=False, index=False, quoting=3)
            logger.info(f"Sorted BED written to: {sorted_bed_file}")

            bigbed_file = self.converter.output_dir / self.converter.genome / f"{self.converter.genome}_sqanti3.bb"
            as_path = os.path.join(self.converter.temp_dir, 'sqanti3_schema.as')
            self._generate_autosql_schema(self.converter.extra_cols, as_path)
            num_extra = len(self.converter.extra_cols)
            cmd = ['bedToBigBed', '-tab', f'-as={as_path}', f'-type=bed12+{num_extra}', '-extraIndex=name', sorted_bed_file, chrom_sizes_file, str(bigbed_file)]
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"BigBed file created: {bigbed_file}")
            return bigbed_file
        except Exception as e:
            logger.error(f"Error creating bigBed file: {e}")
            if result is not None and hasattr(result, 'stderr'):
                logger.error(f"bedToBigBed stderr: {result.stderr}")
            return None
