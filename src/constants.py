"""
Centralized constants for SQANTI-browser: filter definitions and color palettes.
"""

# ---------------------------------------------------------------------------
# BED12 base column names (standard BED columns before any extra fields)
# ---------------------------------------------------------------------------
BED12_BASE_COLS = (
    'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
    'thickStart', 'thickEnd', 'itemRgb', 'blockCount',
    'blockSizes', 'chromStarts',
)

# ---------------------------------------------------------------------------
# Filter definitions for UCSC trackDb and HTML documentation
# Each entry: ('text', col, label) or ('range', col, min, max, label)
# Order matches classification file column order; strand first for trackDb compatibility.
# ---------------------------------------------------------------------------
ORDERED_FILTERS = [
    ('text', 'strand', 'Strand (+ or -)'),
    ('range', 'length', 0, 150000, 'Transcript length (bp)'),
    ('range', 'exons', 0, 400, 'Number of exons'),
    ('text', 'structural_category', 'Structural category (FSM, ISM, NIC, NNC, etc.)'),
    ('text', 'associated_gene', 'Associated gene'),
    ('text', 'associated_transcript', 'Associated transcript'),
    ('range', 'ref_length', 0, 150000, 'Reference transcript length'),
    ('range', 'ref_exons', 0, 400, 'Reference exon count'),
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

# Columns that use dropdown (filterValues) in UCSC; others use text search
DROPDOWN_COLS = frozenset({
    'structural_category', 'subcategory', 'coding', 'FSM_class',
    'all_canonical', 'RTS_stage', 'bite', 'predicted_NMD',
    'within_CAGE_peak', 'polyA_motif_found',
})

# ---------------------------------------------------------------------------
# Default color palettes (RGB tuples)
# ---------------------------------------------------------------------------
DEFAULT_STANDARD_PALETTE = {
    "full-splice_match": (107, 174, 214),
    "incomplete-splice_match": (252, 141, 89),
    "novel_in_catalog": (120, 198, 121),
    "novel_not_in_catalog": (214, 47, 75),
    "genic": (150, 150, 150),
    "antisense": (102, 194, 164),
    "fusion": (218, 165, 32),
    "intergenic": (233, 150, 122),
    "genic_intron": (65, 182, 196),
    "NA": (200, 200, 200),
}

# CSS color names for default palettes (used in UCSC track HTML)
DEFAULT_STANDARD_COLOR_NAMES = {
    "full-splice_match": "steelblue",
    "incomplete-splice_match": "coral",
    "novel_in_catalog": "lightgreen",
    "novel_not_in_catalog": "crimson",
    "genic": "darkgray",
    "antisense": "mediumaquamarine",
    "fusion": "goldenrod",
    "intergenic": "darksalmon",
    "genic_intron": "cadetblue",
    "NA": "silver",
}

DEFAULT_HIGHLIGHT_COLOR_NAMES = {
    "full-splice_match": "steelblue",
    "incomplete-splice_match": "peru",
    "novel_in_catalog": "forestgreen",
    "novel_not_in_catalog": "darkred",
    "genic": "dimgray",
    "antisense": "teal",
    "fusion": "darkgoldenrod",
    "intergenic": "sienna",
    "genic_intron": "darkslategray",
}

DEFAULT_HIGHLIGHT_PALETTE = {
    "full-splice_match": (69, 111, 137),
    "incomplete-splice_match": (202, 113, 71),
    "novel_in_catalog": (77, 126, 78),
    "novel_not_in_catalog": (152, 7, 31),
    "genic": (96, 96, 96),
    "antisense": (66, 124, 105),
    "fusion": (139, 106, 21),
    "intergenic": (149, 96, 78),
    "genic_intron": (42, 117, 126),
}

# Abbreviated names for --category-tracks (user input) -> full structural_category values
CATEGORY_ABBREV_TO_FULL = {
    "FSM": "full-splice_match",
    "ISM": "incomplete-splice_match",
    "NIC": "novel_in_catalog",
    "NNC": "novel_not_in_catalog",
    "antisense": "antisense",
    "genic_intron": "genic_intron",
    "genic_genomic": "genic",  # SQANTI3 uses "genic"
    "genic": "genic",
    "intergenic": "intergenic",
    "fusion": "fusion",
}

DEFAULT_VALIDATION_COLORS = {
    "CAGE_peaks": (200, 133, 41),   # Amber/orange - TSS (field standard)
    "polyA_peaks": (17, 76, 145),   # Blue - TTS (field standard)
    "star_sj": (207, 228, 205),     # Green - splice junctions (field standard)
    "reference": (70, 70, 70),      # Gray
}
