"""
Hub and HTML generation: trackDb, genomes.txt, groups.txt, track HTML, Trix index, category bigBeds.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import logging
from typing import Any, Optional

import pandas as pd

from src.constants import (
    BED12_BASE_COLS,
    DROPDOWN_COLS,
    ORDERED_FILTERS,
    DEFAULT_STANDARD_PALETTE,
    DEFAULT_STANDARD_COLOR_NAMES,
    DEFAULT_HIGHLIGHT_PALETTE,
    DEFAULT_HIGHLIGHT_COLOR_NAMES,
)
from src.utils import normalize_trix_token

logger = logging.getLogger(__name__)


class HubGenerator:
    """Creates UCSC hub files, track HTML, Trix index, and category-specific bigBeds."""

    def __init__(self, converter: Any) -> None:
        self.converter = converter

    def _get_relative_path(self, filename, subdir=None):
        """Generate a relative path for a file in the hub structure."""
        if subdir:
            return f"{subdir}/{filename}"
        return filename

    def _get_rgb_color(self, structural_category):
        """Get RGB color for structural category - uses custom palette if available"""
        palette = (getattr(self.converter, 'final_standard_palette', None) and self.converter.final_standard_palette) or DEFAULT_STANDARD_PALETTE
        cat_key = structural_category.lower()
        rgb_tuple = palette.get(cat_key, (255, 255, 255))
        return f"{rgb_tuple[0]},{rgb_tuple[1]},{rgb_tuple[2]}"

    def _get_category_hex_color(self, structural_category):
        """Get hex color for structural category - uses custom palette if available"""
        palette = (getattr(self.converter, 'final_standard_palette', None) and self.converter.final_standard_palette) or DEFAULT_STANDARD_PALETTE
        cat_key = structural_category.lower().replace(' ', '_')
        rgb_tuple = palette.get(cat_key, (107, 174, 214))
        return f"#{rgb_tuple[0]:02X}{rgb_tuple[1]:02X}{rgb_tuple[2]:02X}"

    def _generate_color_legend_html(self):
        """Generate HTML color legend with color names, hex, and RGB for each structural category."""
        categories = [
            ("full-splice_match", "Full-splice match (FSM)"),
            ("incomplete-splice_match", "Incomplete-splice match (ISM)"),
            ("novel_in_catalog", "Novel in catalog (NIC)"),
            ("novel_not_in_catalog", "Novel not in catalog (NNC)"),
            ("genic", "Genic"),
            ("antisense", "Antisense"),
            ("fusion", "Fusion"),
            ("intergenic", "Intergenic"),
            ("genic_intron", "Genic intron"),
        ]
        standard = getattr(self.converter, 'final_standard_palette', None) or DEFAULT_STANDARD_PALETTE
        highlight = getattr(self.converter, 'final_highlight_palette', None) or DEFAULT_HIGHLIGHT_PALETTE

        def rgb_to_hex(rgb):
            return f"#{rgb[0]:02X}{rgb[1]:02X}{rgb[2]:02X}"

        def get_color_name(cat_key, rgb, palette, name_map):
            default_rgb = palette.get(cat_key)
            if default_rgb and rgb == default_rgb:
                return name_map.get(cat_key, '-')
            return '-'

        html = ['<h2>Color Legend</h2>', '<h3>Standard Colors</h3>', '<table>',
                '<tr><th>Category</th><th>Color Name</th><th>Hex</th><th>RGB</th></tr>']
        for cat_key, cat_name in categories:
            rgb = standard.get(cat_key, (200, 200, 200))
            hex_color = rgb_to_hex(rgb)
            color_name = get_color_name(cat_key, tuple(rgb), DEFAULT_STANDARD_PALETTE, DEFAULT_STANDARD_COLOR_NAMES)
            html.append(f'<tr><td><span class="color-box" style="background-color: {hex_color};"></span>{cat_name}</td>'
                        f'<td>{color_name}</td><td>{hex_color}</td><td>({rgb[0]}, {rgb[1]}, {rgb[2]})</td></tr>')
        html.append('</table>')
        html.extend(['<h3>Highlight Colors (top FL isoform per group)</h3>', '<table>',
                    '<tr><th>Category</th><th>Color Name</th><th>Hex</th><th>RGB</th></tr>'])
        for cat_key, cat_name in categories:
            rgb = highlight.get(cat_key, (128, 128, 128))
            hex_color = rgb_to_hex(rgb)
            color_name = get_color_name(cat_key, tuple(rgb), DEFAULT_HIGHLIGHT_PALETTE, DEFAULT_HIGHLIGHT_COLOR_NAMES)
            html.append(f'<tr><td><span class="color-box" style="background-color: {hex_color};"></span>{cat_name}</td>'
                        f'<td>{color_name}</td><td>{hex_color}</td><td>({rgb[0]}, {rgb[1]}, {rgb[2]})</td></tr>')
        html.append('</table>')
        return '\n    '.join(html)

    def _generate_color_legend_markdown(self):
        """Generate Markdown color legend using current palette (custom or default)"""
        categories = [
            ("full-splice_match", "Full-splice Match (FSM)"),
            ("incomplete-splice_match", "Incomplete-splice Match (ISM)"),
            ("novel_in_catalog", "Novel In Catalog (NIC)"),
            ("novel_not_in_catalog", "Novel Not In Catalog (NNC)"),
            ("genic", "Genic"),
            ("antisense", "Antisense"),
            ("fusion", "Fusion"),
            ("intergenic", "Intergenic"),
            ("genic_intron", "Genic Intron"),
        ]
        standard = getattr(self.converter, 'final_standard_palette', None) or DEFAULT_STANDARD_PALETTE
        highlight = getattr(self.converter, 'final_highlight_palette', None) or DEFAULT_HIGHLIGHT_PALETTE

        def rgb_to_hex(rgb):
            return f"#{rgb[0]:02X}{rgb[1]:02X}{rgb[2]:02X}"

        lines = ['## 🎨 Color Legend', '', '### Standard Colors', '']
        for cat_key, cat_name in categories:
            rgb = standard.get(cat_key, (200, 200, 200))
            lines.append(f'- **{cat_name}:** {rgb_to_hex(rgb)} (RGB {rgb[0]}, {rgb[1]}, {rgb[2]})')
        lines.extend(['', '### Highlight Colors (top FL isoform per group)', ''])
        for cat_key, cat_name in categories:
            rgb = highlight.get(cat_key, (128, 128, 128))
            lines.append(f'- **{cat_name}:** {rgb_to_hex(rgb)} (RGB {rgb[0]}, {rgb[1]}, {rgb[2]})')
        return '\n'.join(lines)

    def _generate_filter_options_html(self, category_filter=None):
        """Generate HTML documentation for filter options."""
        if not hasattr(self.converter, 'classification_df'):
            return "<p>No classification data available.</p>"
        df = self.converter.classification_df
        if category_filter:
            df = df[df['structural_category'] == category_filter]
        html_parts = []
        html_parts.append('<table>')
        html_parts.append('<tr><th>Filter</th><th>Type</th><th>Values / Range</th></tr>')
        for filter_def in ORDERED_FILTERS:
            if filter_def[0] == 'text':
                _, col, label = filter_def
                if col in df.columns:
                    if col in ['associated_gene', 'associated_transcript', 'min_cov_pos']:
                        unique_count = df[col].nunique()
                        html_parts.append(f'<tr><td><strong>{label}</strong></td><td>🔤 Text</td><td>{unique_count} unique values. Use wildcards like <code>*</code></td></tr>')
                    elif col in DROPDOWN_COLS:
                        value_counts = df[col].value_counts(dropna=False)
                        unique_vals = [f"<code>{val}</code> ({count})" for val, count in value_counts.items() if pd.notna(val) and str(val) not in ('NA', '')]
                        if unique_vals:
                            display_vals = '<br>'.join(unique_vals)
                            html_parts.append(f'<tr><td><strong>{label}</strong></td><td>📋 Dropdown</td><td class="filter-values">{display_vals}</td></tr>')
                    else:
                        value_counts = df[col].value_counts(dropna=False)
                        unique_vals = [f"<code>{val}</code> ({count})" for val, count in value_counts.items() if pd.notna(val) and str(val) not in ('NA', '')]
                        if unique_vals:
                            display_vals = ', '.join(unique_vals[:8]) + f', ... ({len(unique_vals)} total)' if len(unique_vals) > 8 else ', '.join(unique_vals)
                            html_parts.append(f'<tr><td><strong>{label}</strong></td><td>🔤 Text</td><td class="filter-values">{display_vals}</td></tr>')
            else:
                _, col, min_val, max_val, label = filter_def
                if col in df.columns:
                    try:
                        data_min = df[col].dropna().astype(float).min()
                        data_max = df[col].dropna().astype(float).max()
                        html_parts.append(f'<tr><td><strong>{label}</strong></td><td>📊 Slider</td><td>Data range: {data_min:.0f} to {data_max:.0f}</td></tr>')
                    except (ValueError, TypeError):
                        html_parts.append(f'<tr><td><strong>{label}</strong></td><td>📊 Slider</td><td>Range: {min_val} to {max_val}</td></tr>')
        html_parts.append('<tr><td><strong>Number of exons (from BED)</strong></td><td>📊 Slider</td><td>Range: 0 to 400</td></tr>')
        html_parts.append('</table>')
        return '\n    '.join(html_parts)

    def _generate_trix_index(self, enhanced_bed_file, genome_dir):
        """Generate Trix (.ix/.ixx) text index with rich descriptions for fast search"""
        if not shutil.which('ixIxx'):
            logger.warning("ixIxx not found; skipping Trix index generation")
            return False
        logger.info("Generating Trix index...")
        try:
            trix_input = os.path.join(self.converter.temp_dir, 'trix_input.txt')
            if not hasattr(self.converter, 'classification_df'):
                logger.error("Classification data not loaded")
                return False
            classification_data = self.converter.classification_df.set_index('name').to_dict('index')
            with open(enhanced_bed_file, 'r') as bed_in, open(trix_input, 'w') as t_out:
                for line in bed_in:
                    parts = line.rstrip('\t\n').split('\t')
                    if len(parts) >= 12:
                        tid = parts[3]
                        d = classification_data.get(tid, {})
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
                        desc = f"{tid} " + ' '.join([f"{k} {v}" for k, v in fields.items() if v not in (None, '', 'NA', 'nan')])
                        synonym_parts = []
                        synonym_seen = set()

                        def add_synonym(token):
                            if token and token not in synonym_seen:
                                synonym_parts.append(token)
                                synonym_seen.add(token)

                        for k, v in fields.items():
                            if v not in (None, '', 'NA', 'nan'):
                                v_str = str(v)
                                add_synonym(v_str)
                                add_synonym(f"{k}_{v_str}")
                                normalized_key = normalize_trix_token(k)
                                normalized_value = normalize_trix_token(v)
                                if normalized_value:
                                    add_synonym(normalized_value)
                                if normalized_key and normalized_value:
                                    add_synonym(f"{normalized_key}_{normalized_value}")
                        synonyms = ' '.join(synonym_parts)
                        t_out.write(f"{tid}\t{desc}\t{synonyms}\n")

            ix_path = os.path.join(genome_dir, 'trix.ix')
            ixx_path = os.path.join(genome_dir, 'trix.ixx')
            cmd = ['ixIxx', '-maxWordLength=64', trix_input, ix_path, ixx_path]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"Trix index generated: {ix_path}, {ixx_path}")
            return True
        except Exception as e:
            logger.warning(f"Failed to generate Trix index: {e}")
            import traceback
            logger.warning(traceback.format_exc())
            return False

    def create_category_bigbeds(self, main_bed_file: str) -> bool:
        """Create separate bigBed files for each structural category"""
        logger.info("Creating category-specific bigBed files...")
        try:
            as_path = os.path.join(self.converter.temp_dir, 'sqanti3_schema.as')
            num_extra = len(self.converter.extra_cols)
            chrom_sizes_file = self.converter.chrom_sizes_file or os.path.join(self.converter.temp_dir, 'chrom.sizes')
            bed_cols = list(BED12_BASE_COLS) + self.converter.extra_cols
            df = pd.read_csv(main_bed_file, sep='\t', names=bed_cols, dtype={'chrom': 'string'})
            if 'structural_category' not in df.columns:
                logger.warning("structural_category column missing, cannot create category tracks")
                return False
            categories = df['structural_category'].unique()
            self.converter.category_bigbeds = {}
            bed_processor = getattr(self.converter, '_bed_processor', None)
            if bed_processor is None:
                logger.error("BedProcessor not attached to converter")
                return False
            for cat in categories:
                if pd.isna(cat) or cat == 'NA':
                    continue
                cat_str = str(cat)
                safe_cat_file = cat_str.replace(' ', '_').replace('/', '_')
                sub_sorted = os.path.join(self.converter.temp_dir, f"sqanti3_{safe_cat_file}.sorted.bed")
                sub_bb = self.converter.output_dir / self.converter.genome / f"{self.converter.genome}_sqanti3_{safe_cat_file}.bb"
                cat_df = df[df['structural_category'] == cat].copy()
                cat_df_sorted = bed_processor.sort_bed_for_visualization(cat_df, self.converter.sort_by)
                cat_df_sorted.to_csv(sub_sorted, sep='\t', header=False, index=False, quoting=3)
                cmd = ['bedToBigBed', '-tab', f'-as={as_path}', f'-type=bed12+{num_extra}', '-extraIndex=name', sub_sorted, chrom_sizes_file, str(sub_bb)]
                result = subprocess.run(cmd, check=False, capture_output=True, text=True)
                if result.returncode != 0:
                    logger.error(f"bedToBigBed failed for category {cat}")
                    logger.error(f"Command: {' '.join(cmd)}")
                    logger.error(f"Stderr: {result.stderr}")
                    raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
                self.converter.category_bigbeds[cat_str] = sub_bb.name
                logger.info(f"Created category track: {sub_bb}")
            return True
        except Exception as e:
            logger.error(f"Error creating category bigBeds: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False

    def _create_hub_txt(self, hub_name: str) -> None:
        """Create hub.txt."""
        c = self.converter
        with open(c.output_dir / "hub.txt", 'w', newline='\n') as f:
            f.write(f"hub {hub_name}\n")
            f.write(f"shortLabel {hub_name}\n")
            f.write(f"longLabel SQANTI3 Transcriptome Analysis for {c.genome}\n")
            f.write(f"genomesFile {self._get_relative_path('genomes.txt')}\n")
            f.write(f"email sqanti3_user@users.noreply.github.com\n")
            f.write(f"descriptionUrl {self._get_relative_path('README.md')}\n")

    def _create_genomes_txt(self, hub_name: str) -> None:
        """Create genomes.txt (and optional 2bit genome description HTML)."""
        c = self.converter
        genomes_file = c.output_dir / "genomes.txt"
        with open(genomes_file, 'w', newline='\n') as f:
            f.write(f"genome {c.genome}\n")
            f.write(f"trackDb {self._get_relative_path('trackDb.txt', subdir=c.genome)}\n")
            f.write(f"groups {self._get_relative_path('groups.txt', subdir=c.genome)}\n")
            if c.two_bit_file:
                genome_2bit_name = f"{c.genome}.2bit"
                shutil.copy(c.two_bit_file, c.output_dir / c.genome / genome_2bit_name)
                f.write(f"twoBitPath {self._get_relative_path(genome_2bit_name, subdir=c.genome)}\n")
                f.write(f"organism user-defined\n")
                f.write(f"description user-defined\n")
                f.write(f"scientificName user-defined\n")
                f.write(f"orderKey 4800\n")
                genome_html_name = f"{c.genome}_genome.html"
                with open(c.output_dir / c.genome / genome_html_name, 'w') as gf:
                    gf.write(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{c.genome} Description</title>
    <style>
        body {{ font-family: sans-serif; margin: 20px; }}
        h1 {{ color: #006699; }}
        .info-box {{ background: #f0f0f0; padding: 15px; border-radius: 5px; }}
    </style>
</head>
<body>
    <h1>User-Defined Genome: {c.genome}</h1>
    <div class="info-box">
        <p><strong>Organism:</strong> User-defined</p>
        <p><strong>Description:</strong> Custom reference genome provided via 2bit file.</p>
        <p><strong>Source File:</strong> {os.path.basename(c.two_bit_file)}</p>
    </div>
    <h2>Chromosomes</h2>
    <ul>""")
                    try:
                        with open(os.path.join(c.temp_dir, 'chrom.sizes'), 'r') as cf:
                            for line in cf:
                                parts = line.strip().split('\t')
                                if len(parts) >= 2:
                                    gf.write(f"<li><strong>{parts[0]}:</strong> {parts[1]} bp</li>\n")
                    except Exception:
                        gf.write("<li>Could not read chromosome sizes.</li>")
                    gf.write("""    </ul>
</body>
</html>""")
                f.write(f"htmlPath {self._get_relative_path(genome_html_name, subdir=c.genome)}\n")
                try:
                    with open(os.path.join(c.temp_dir, 'chrom.sizes'), 'r') as cf:
                        first_line = cf.readline().strip()
                        if first_line:
                            chrom = first_line.split('\t')[0]
                            f.write(f"defaultPos {chrom}:1-5000\n")
                except Exception as e:
                    logger.warning(f"Could not determine default position for genomes.txt: {e}")

    def _create_groups_txt(self) -> None:
        """Create groups.txt in genome directory."""
        c = self.converter
        genome_dir = c.output_dir / c.genome
        genome_dir.mkdir(exist_ok=True)
        with open(genome_dir / "groups.txt", 'w', newline='\n') as gf:
            gf.write("name transcripts\n")
            gf.write("label Transcripts\n")
            gf.write("priority 1\n")
            gf.write("defaultIsClosed 0\n\n")
            gf.write("name validation\n")
            gf.write("label Validation Data\n")
            gf.write("priority 2\n")
            gf.write("defaultIsClosed 0\n")

    def _create_transcripts_track_html(self, filter_options_html: str, color_legend_html: str) -> str:
        """Create main SQANTI3 transcripts track HTML. Returns filename."""
        c = self.converter
        transcripts_html_name = f"{c.genome}_sqanti3_track.html"
        with open(c.output_dir / c.genome / transcripts_html_name, 'w') as tf:
            tf.write(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{c.genome} SQANTI3 Transcripts</title>
    <style>
        html, body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif; line-height: 1.6; margin: 12px !important; background-color: transparent !important; color: #212529; text-align: left !important; }}
        .sqanti-wrap {{ max-width: 900px; margin: 0; background: transparent; padding: 0; }}
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
<div class="sqanti-wrap">
    <h1>SQANTI3 Transcripts ({c.genome})</h1>
    <p>This track displays <strong>{len(c.classification_df)}</strong> SQANTI3 transcript models. Colors indicate structural category.</p>
    <h2>How to Filter</h2>
    <div class="tip">
        <strong>Tip:</strong> Right-click on the track and select "Configure" to access all filters.
        Use <code>*</code> as a wildcard in text filters (e.g., <code>FSM*</code> matches FSMA, FSMB, FSMC, FSMD).
    </div>
    <h2>Available Filters</h2>
    <p>All filters available for this track:</p>
{filter_options_html}
    {color_legend_html}
</div>
</body>
</html>""")
        return transcripts_html_name

    def _create_validation_track_html_files(self) -> dict:
        """Create HTML for validation tracks (star, cage, polya, ref). Returns dict of track_name -> html_filename."""
        c = self.converter
        result = {}
        if c.star_bigbed and os.path.exists(c.star_bigbed):
            name = f"{c.genome}_star_sj_track.html"
            with open(c.output_dir / c.genome / name, 'w') as f:
                f.write(f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8"><title>{c.genome} STAR Splice Junctions</title></head>
<body><h1>STAR Splice Junctions ({c.genome})</h1>
<p>Junctions converted from STAR <code>SJ.out.tab</code> (bed6/bigBed).</p></body>
</html>""")
            result['star'] = name
        if c.cage_bigbed and os.path.exists(c.cage_bigbed):
            name = f"{c.genome}_cage_peaks_track.html"
            with open(c.output_dir / c.genome / name, 'w') as f:
                f.write(f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8"><title>{c.genome} CAGE Peaks</title></head>
<body><h1>CAGE Peaks ({c.genome})</h1>
<p>Cap Analysis of Gene Expression (CAGE) peaks, representing transcription start sites (TSS).</p></body>
</html>""")
            result['cage'] = name
        if c.polya_bigbed and os.path.exists(c.polya_bigbed):
            name = f"{c.genome}_polya_peaks_track.html"
            with open(c.output_dir / c.genome / name, 'w') as f:
                f.write(f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8"><title>{c.genome} PolyA Peaks</title></head>
<body><h1>PolyA Peaks ({c.genome})</h1>
<p>PolyA site peaks, representing transcription termination sites (TTS).</p></body>
</html>""")
            result['polya'] = name
        if c.ref_bigbed and os.path.exists(c.ref_bigbed):
            name = f"{c.genome}_reference_track.html"
            with open(c.output_dir / c.genome / name, 'w') as f:
                f.write(f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8"><title>{c.genome} Reference Annotation</title></head>
<body><h1>Reference Annotation ({c.genome})</h1>
<p>Reference genome annotation used for SQANTI3 analysis.</p>
<p><strong>Source:</strong> {os.path.basename(c.ref_gtf)}</p></body>
</html>""")
            result['ref'] = name
        return result

    def _write_trackdb(self, transcripts_html_name: str, validation_html: dict, num_extra: int) -> None:
        """Write trackDb.txt with main track, validation tracks, and category tracks."""
        c = self.converter
        genome_dir = c.output_dir / c.genome

        def sanitize(col):
            return col.replace('.', '_').replace(' ', '_').replace('-', '_').replace('/', '_')

        existing_fields = set(sanitize(col) for col in c.extra_cols)

        def get_unique_values(col):
            if col in c.classification_df.columns:
                vals = c.classification_df[col].dropna().astype(str).unique()
                return sorted([v for v in vals if v and v != 'NA' and v != 'nan'])
            return []

        trackdb_file = genome_dir / "trackDb.txt"
        star_html_name = validation_html.get('star')
        cage_html_name = validation_html.get('cage')
        polya_html_name = validation_html.get('polya')
        ref_html_name = validation_html.get('ref')
        with open(trackdb_file, 'w', newline='\n') as f:
            f.write(f"track {c.genome}_sqanti3\n")
            f.write(f"bigDataUrl {self._get_relative_path(f'{c.genome}_sqanti3.bb')}\n")
            f.write(f"shortLabel SQANTI3 Transcripts\n")
            f.write(f"longLabel SQANTI3 Transcriptome Analysis Results\n")
            f.write(f"type bigBed 12 + {num_extra}\n")
            for filter_def in ORDERED_FILTERS:
                if filter_def[0] == 'text':
                    _, col, label = filter_def
                    if col in existing_fields:
                        if col in DROPDOWN_COLS:
                            unique_vals = get_unique_values(col)
                            if unique_vals:
                                f.write(f"filterValues.{col} {','.join(unique_vals)}\n")
                                f.write(f"filterLabel.{col} {label}\n")
                            else:
                                f.write(f"filterText.{col} *\n")
                                f.write(f"filterType.{col} wildcard\n")
                                f.write(f"filterLabel.{col} {label}\n")
                        else:
                            f.write(f"filterText.{col} *\n")
                            f.write(f"filterType.{col} wildcard\n")
                            f.write(f"filterLabel.{col} {label}\n")
                else:
                    _, col, min_val, max_val, label = filter_def
                    if col in existing_fields:
                        f.write(f"filter.{col} {min_val}:{max_val}\n")
                        f.write(f"filterByRange.{col} on\n")
                        f.write(f"filterLimits.{col} {min_val}:{max_val}\n")
                        f.write(f"filterLabel.{col} {label}\n")
            f.write(f"filter.blockCount 0:400\n")
            f.write(f"filterByRange.blockCount on\n")
            f.write(f"filterLimits.blockCount 0:400\n")
            f.write(f"filterLabel.blockCount Number of exons (from BED)\n")
            f.write(f"visibility pack\n")
            f.write(f"group transcripts\n")
            f.write(f"itemRgb on\n")
            f.write(f"priority 1\n")
            f.write(f"html {self._get_relative_path(transcripts_html_name)}\n")
            f.write(f"searchIndex name\n")
            trix_ix = genome_dir / 'trix.ix'
            if os.path.exists(trix_ix):
                f.write(f"searchTrix trix.ix\n")

            if c.star_bigbed and os.path.exists(c.star_bigbed):
                f.write("\ntrack " + f"{c.genome}_star_sj\n")
                f.write(f"bigDataUrl {self._get_relative_path(os.path.basename(c.star_bigbed))}\n")
                f.write(f"shortLabel STAR Junctions\n")
                f.write(f"longLabel STAR splice junctions (SJ.out.tab)\n")
                f.write(f"type bigBed 9\n")
                f.write(f"itemRgb on\nvisibility full\ngroup validation\npriority 2\n")
                if star_html_name:
                    f.write(f"html {self._get_relative_path(star_html_name)}\n")
            if c.cage_bigbed and os.path.exists(c.cage_bigbed):
                f.write("\ntrack " + f"{c.genome}_cage_peaks\n")
                f.write(f"bigDataUrl {self._get_relative_path(os.path.basename(c.cage_bigbed))}\n")
                f.write(f"shortLabel CAGE Peaks\nlongLabel CAGE Peaks (TSS Validation)\n")
                f.write(f"type bigBed 9\nitemRgb on\nvisibility dense\ngroup validation\npriority 3\n")
                if cage_html_name:
                    f.write(f"html {self._get_relative_path(cage_html_name)}\n")
            if c.polya_bigbed and os.path.exists(c.polya_bigbed):
                f.write("\ntrack " + f"{c.genome}_polya_peaks\n")
                f.write(f"bigDataUrl {self._get_relative_path(os.path.basename(c.polya_bigbed))}\n")
                f.write(f"shortLabel PolyA Peaks\nlongLabel PolyA Site Peaks (TTS Validation)\n")
                f.write(f"type bigBed 9\nitemRgb on\nvisibility dense\ngroup validation\npriority 4\n")
                if polya_html_name:
                    f.write(f"html {self._get_relative_path(polya_html_name)}\n")
            if c.ref_bigbed and os.path.exists(c.ref_bigbed):
                f.write("\ntrack " + f"{c.genome}_reference\n")
                f.write(f"bigDataUrl {self._get_relative_path(os.path.basename(c.ref_bigbed))}\n")
                f.write(f"shortLabel Reference\nlongLabel Reference Genome Annotation\n")
                f.write(f"type bigBed 12\nvisibility pack\ngroup validation\npriority 5\ncolor 70,70,70\n")
                if ref_html_name:
                    f.write(f"html {self._get_relative_path(ref_html_name)}\n")

            if hasattr(c, 'category_bigbeds') and c.category_bigbeds:
                for cat, bb_filename in c.category_bigbeds.items():
                    safe_cat = cat.replace(' ', '_').replace('/', '_')
                    cat_html_name = f"{c.genome}_sqanti3_{safe_cat}.html"
                    count = len(c.classification_df[c.classification_df['structural_category'] == cat]) if hasattr(c, 'classification_df') else 0
                    cat_filter_options = self._generate_filter_options_html(category_filter=cat)
                    cat_html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"><title>{c.genome} SQANTI3 {cat} Transcripts</title>
<style>
html, body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif; margin: 12px !important; background-color: transparent !important; color: #212529; }}
.sqanti-wrap {{ max-width: 900px; }}
h1 {{ color: #333; border-bottom: 2px solid {self._get_category_hex_color(cat)}; }}
table {{ border-collapse: collapse; width: 100%; }}
th, td {{ border: 1px solid #dee2e6; padding: 10px; }}
.category-badge {{ display: inline-block; background-color: {self._get_category_hex_color(cat)}; color: white; padding: 5px 12px; border-radius: 15px; font-weight: bold; }}
</style>
</head>
<body><div class="sqanti-wrap">
<h1>SQANTI3 Transcripts: <span class="category-badge">{cat}</span></h1>
<p>This track contains <strong>{count}</strong> transcripts classified as <strong>{cat}</strong> for <strong>{c.genome}</strong>.</p>
<h2>Available Filters</h2>
{cat_filter_options}
</div></body>
</html>"""
                    try:
                        with open(c.output_dir / c.genome / cat_html_name, 'w') as ch:
                            ch.write(cat_html_content)
                    except Exception as e:
                        logger.error(f"Error creating HTML for category {cat}: {e}")
                    f.write("\ntrack " + f"{c.genome}_sqanti3_{safe_cat}\n")
                    f.write(f"bigDataUrl {self._get_relative_path(os.path.basename(bb_filename))}\n")
                    f.write(f"shortLabel SQANTI3 {cat}\n")
                    f.write(f"longLabel SQANTI3 {cat} transcripts\n")
                    f.write(f"type bigBed 12 + {num_extra}\n")
                    f.write(f"visibility hide\ngroup transcripts\nitemRgb on\npriority 3\n")
                    f.write(f"html {self._get_relative_path(cat_html_name)}\n")
                    cat_df = c.classification_df[c.classification_df['structural_category'] == cat]

                    def get_cat_unique_values(col):
                        if col in cat_df.columns:
                            vals = cat_df[col].dropna().astype(str).unique()
                            return sorted([v for v in vals if v and v != 'NA' and v != 'nan'])
                        return []

                    for filter_def in ORDERED_FILTERS:
                        if filter_def[0] == 'text':
                            _, col, label = filter_def
                            if col in existing_fields:
                                if col in DROPDOWN_COLS:
                                    unique_vals = get_cat_unique_values(col)
                                    if unique_vals:
                                        f.write(f"filterValues.{col} {','.join(unique_vals)}\n")
                                        f.write(f"filterLabel.{col} {label}\n")
                                    else:
                                        f.write(f"filterText.{col} *\nfilterType.{col} wildcard\nfilterLabel.{col} {label}\n")
                                else:
                                    f.write(f"filterText.{col} *\nfilterType.{col} wildcard\nfilterLabel.{col} {label}\n")
                        else:
                            _, col, min_val, max_val, label = filter_def
                            if col in existing_fields:
                                f.write(f"filter.{col} {min_val}:{max_val}\n")
                                f.write(f"filterByRange.{col} on\nfilterLimits.{col} {min_val}:{max_val}\n")
                                f.write(f"filterLabel.{col} {label}\n")
                    f.write(f"filter.blockCount 0:400\nfilterByRange.blockCount on\n")
                    f.write(f"filterLimits.blockCount 0:400\nfilterLabel.blockCount Number of exons (from BED)\n")

    def _create_readme(self, hub_name: str, num_extra: int) -> None:
        """Create README.md for the hub."""
        c = self.converter
        color_legend_md = self._generate_color_legend_markdown()
        with open(c.output_dir / "README.md", 'w') as f_md:
            f_md.write(f"""# {hub_name}

This hub displays SQANTI3 transcriptome analysis results for the {c.genome} genome assembly.

## 🚀 Usage Instructions

1. Upload all files in this directory to a web-accessible location (e.g., GitHub).
2. In the UCSC Genome Browser, go to **My Data → Track Hubs**.
3. Enter the URL to your `hub.txt` file.
4. Select the appropriate genome assembly ({c.genome}).
5. The SQANTI3 tracks will appear in your track list.

## 🔍 Advanced Filtering

This hub uses the bigBed 12+{num_extra} format with native UCSC filters. Right-click on the track and select "Configure" or "Filter" to access controls.

## 🔎 Trix Search

Use the search box to find isoforms by attribute. Search terms use underscores (e.g., `structural_category_full_splice_match`).

{color_legend_md}
""")

    def create_hub_files(self, bigbed_file: Any) -> bool:
        """Create UCSC Genome Browser hub files."""
        logger.info("Creating hub files...")
        c = self.converter
        output_dir_name = c.output_dir.name
        hub_name = f"{output_dir_name}_{c.genome}_SQANTI3_Hub"
        num_extra = len(c.extra_cols) if hasattr(c, 'extra_cols') else 0

        self._create_hub_txt(hub_name)
        self._create_groups_txt()  # Creates genome_dir (needed before genomes.txt 2bit copy)
        self._create_genomes_txt(hub_name)
        filter_options_html = self._generate_filter_options_html()
        color_legend_html = self._generate_color_legend_html()
        transcripts_html_name = self._create_transcripts_track_html(filter_options_html, color_legend_html)
        validation_html = self._create_validation_track_html_files()
        self._write_trackdb(transcripts_html_name, validation_html, num_extra)
        self._create_readme(hub_name, num_extra)

        logger.info("Hub files created successfully")
        return True
