import pandas as pd
import argparse
import os
import json
from pathlib import Path
import sys
import logging
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Default color palettes (RGB tuples)
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
}

def rgb_to_hex(rgb):
    """Convert RGB tuple to hex string."""
    return f"#{rgb[0]:02X}{rgb[1]:02X}{rgb[2]:02X}"


def normalize_trix_token(value):
    """Normalize strings for Trix search tokens (lowercase, underscores)."""
    if value is None:
        return ""
    token = str(value).strip().lower()
    if not token:
        return ""
    if token == '+':
        return 'plus'
    if token == '-':
        return 'minus'
    token = re.sub(r'[^\w]+', '_', token)
    token = re.sub(r'_+', '_', token)
    token = token.strip('_')
    return token


def get_category_color(category, palette=None):
    """Get hex color for a category, using custom palette if provided."""
    if palette is None:
        palette = DEFAULT_STANDARD_PALETTE
    rgb = palette.get(category, (200, 200, 200))
    return rgb_to_hex(rgb)


def generate_category_svg(category, palette=None, include_reference=True):
    """Generate SVG for a structural category with custom color."""
    color = get_category_color(category, palette)
    
    # SVG templates with placeholders for color
    if include_reference:
        # Full SVGs with reference (for individual category reports)
        templates = {
            "full-splice_match": f'''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="190" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="80" y="8" width="40" height="8" fill="black" rx="2"/><rect x="150" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><line x1="10" y1="52" x2="190" y2="52" stroke="{color}" stroke-width="2" /><rect x="10" y="48" width="40" height="8" fill="{color}" rx="2"/><rect x="80" y="48" width="40" height="8" fill="{color}" rx="2"/><rect x="150" y="48" width="40" height="8" fill="{color}" rx="2"/><text x="220" y="55" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (FSM)</text></svg>''',
            "incomplete-splice_match": f'''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="190" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="80" y="8" width="40" height="8" fill="black" rx="2"/><rect x="150" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><line x1="80" y1="52" x2="190" y2="52" stroke="{color}" stroke-width="2" /><rect x="80" y="48" width="40" height="8" fill="{color}" rx="2"/><rect x="150" y="48" width="40" height="8" fill="{color}" rx="2"/><text x="220" y="55" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (ISM)</text></svg>''',
            "novel_in_catalog": f'''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="190" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="80" y="8" width="40" height="8" fill="black" rx="2"/><rect x="150" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><line x1="10" y1="52" x2="190" y2="52" stroke="{color}" stroke-width="2" /><rect x="10" y="48" width="40" height="8" fill="{color}" rx="2"/><rect x="150" y="48" width="40" height="8" fill="{color}" rx="2"/><text x="220" y="55" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (NIC)</text></svg>''',
            "novel_not_in_catalog": f'''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="190" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="80" y="8" width="40" height="8" fill="black" rx="2"/><rect x="150" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><line x1="10" y1="52" x2="190" y2="52" stroke="{color}" stroke-width="2" /><rect x="10" y="48" width="40" height="8" fill="{color}" rx="2"/><rect x="80" y="48" width="60" height="8" fill="{color}" rx="2"/><rect x="150" y="48" width="40" height="8" fill="{color}" rx="2"/><text x="130" y="40" font-size="10" fill="{color}" text-anchor="middle">New Site</text><text x="220" y="55" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (NNC)</text></svg>''',
            "genic_intron": f'''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="190" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="150" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><rect x="80" y="48" width="40" height="8" fill="{color}" rx="2"/><text x="220" y="55" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (Intron)</text></svg>''',
            "genic": f'''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="120" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="80" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><rect x="40" y="48" width="60" height="8" fill="{color}" rx="2"/><text x="70" y="70" font-size="10" fill="{color}" text-anchor="middle">Overlap</text><text x="220" y="55" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (Genic)</text></svg>''',
            "antisense": f'''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="120" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="80" y="8" width="40" height="8" fill="black" rx="2"/><text x="130" y="15" font-size="12">→</text><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference (+)</text><line x1="10" y1="52" x2="120" y2="52" stroke="{color}" stroke-width="2" /><rect x="10" y="48" width="40" height="8" fill="{color}" rx="2"/><rect x="80" y="48" width="40" height="8" fill="{color}" rx="2"/><text x="130" y="55" font-size="12" fill="{color}">←</text><text x="220" y="55" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (-)</text></svg>''',
            "intergenic": f'''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><rect x="140" y="48" width="40" height="8" fill="{color}" rx="2"/><text x="220" y="55" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (Inter)</text></svg>''',
            "fusion": f'''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><text x="55" y="15" font-size="10">Gene A</text><rect x="140" y="8" width="40" height="8" fill="black" rx="2"/><text x="185" y="15" font-size="10">Gene B</text><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><line x1="10" y1="52" x2="180" y2="52" stroke="{color}" stroke-width="2" /><rect x="10" y="48" width="40" height="8" fill="{color}" rx="2"/><rect x="140" y="48" width="40" height="8" fill="{color}" rx="2"/><text x="220" y="55" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (Fusion)</text></svg>'''
        }
    else:
        # Overview SVGs without reference (for complete transcriptome report)
        templates = {
            "full-splice_match": f'''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="30" x2="190" y2="30" stroke="{color}" stroke-width="2" /><rect x="10" y="26" width="40" height="8" fill="{color}" rx="2"/><rect x="80" y="26" width="40" height="8" fill="{color}" rx="2"/><rect x="150" y="26" width="40" height="8" fill="{color}" rx="2"/><text x="220" y="33" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (FSM)</text></svg>''',
            "incomplete-splice_match": f'''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><line x1="80" y1="30" x2="190" y2="30" stroke="{color}" stroke-width="2" /><rect x="80" y="26" width="40" height="8" fill="{color}" rx="2"/><rect x="150" y="26" width="40" height="8" fill="{color}" rx="2"/><text x="220" y="33" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (ISM)</text></svg>''',
            "novel_in_catalog": f'''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="30" x2="190" y2="30" stroke="{color}" stroke-width="2" /><rect x="10" y="26" width="40" height="8" fill="{color}" rx="2"/><rect x="150" y="26" width="40" height="8" fill="{color}" rx="2"/><text x="220" y="33" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (NIC)</text></svg>''',
            "novel_not_in_catalog": f'''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="30" x2="190" y2="30" stroke="{color}" stroke-width="2" /><rect x="10" y="26" width="40" height="8" fill="{color}" rx="2"/><rect x="80" y="26" width="60" height="8" fill="{color}" rx="2"/><rect x="150" y="26" width="40" height="8" fill="{color}" rx="2"/><text x="130" y="18" font-size="10" fill="{color}" text-anchor="middle">New Site</text><text x="220" y="33" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (NNC)</text></svg>''',
            "genic_intron": f'''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><rect x="80" y="26" width="40" height="8" fill="{color}" rx="2"/><text x="220" y="33" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (Intron)</text></svg>''',
            "genic": f'''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><rect x="40" y="26" width="60" height="8" fill="{color}" rx="2"/><text x="70" y="48" font-size="10" fill="{color}" text-anchor="middle">Overlap</text><text x="220" y="33" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (Genic)</text></svg>''',
            "antisense": f'''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="30" x2="120" y2="30" stroke="{color}" stroke-width="2" /><rect x="10" y="26" width="40" height="8" fill="{color}" rx="2"/><rect x="80" y="26" width="40" height="8" fill="{color}" rx="2"/><text x="130" y="33" font-size="12" fill="{color}">←</text><text x="220" y="33" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (-)</text></svg>''',
            "intergenic": f'''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><rect x="140" y="26" width="40" height="8" fill="{color}" rx="2"/><text x="220" y="33" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (Inter)</text></svg>''',
            "fusion": f'''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="30" x2="180" y2="30" stroke="{color}" stroke-width="2" /><rect x="10" y="26" width="40" height="8" fill="{color}" rx="2"/><rect x="140" y="26" width="40" height="8" fill="{color}" rx="2"/><text x="220" y="33" font-family="sans-serif" font-size="12" fill="{color}" font-weight="bold">Isoform (Fusion)</text></svg>'''
        }
    
    return templates.get(category, "")

# HTML Template with DataTables
# We use CDNs for jQuery and DataTables. 
# Includes column-specific filtering in the footer and export buttons.
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{report_title}</title>
    
    <!-- DataTables CSS -->
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/buttons/2.4.2/css/buttons.dataTables.min.css">
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/fixedheader/3.4.0/css/fixedHeader.dataTables.min.css">
    
    <style>
        body {{ 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
            margin: 20px; 
            background-color: #f8f9fa;
        }}
        h2 {{ color: #2c3e50; }}
        .container {{ 
            background-color: white; 
            padding: 20px; 
            border-radius: 8px; 
            box-shadow: 0 2px 4px rgba(0,0,0,0.1); 
        }}
        .info-box {{
            background-color: #f1f3f5;
            padding: 20px;
            border-radius: 6px;
            margin-bottom: 20px;
            border-left: 5px solid #6c757d;
        }}
        .info-box p {{
            margin-top: 5px;
            margin-bottom: 5px;
        }}
        table.dataTable thead th {{
            background-color: #e9ecef;
            color: #495057;
            white-space: nowrap;
        }}
        /* Table cells */
        table.dataTable tbody td {{
            white-space: nowrap;
            font-size: 0.9em;
            padding: 4px 10px;
        }}
        /* Filter input styling */
        tfoot input, thead select {{
            width: 100%;
            padding: 3px;
            box-sizing: border-box;
            border: 1px solid #ced4da;
            border-radius: 4px;
            font-size: 0.8em;
        }}
        .dataTables_wrapper {{ margin-top: 20px; }}
        /* Row selection styling */
        table.dataTable tbody tr.selected {{
            background-color: #b3d4fc !important;
        }}
        table.dataTable tbody tr:hover {{
            cursor: pointer;
        }}
    </style>

    <!-- jQuery and DataTables JS -->
    <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/2.4.2/js/dataTables.buttons.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.html5.min.js"></script>
    <script src="https://cdn.datatables.net/fixedheader/3.4.0/js/dataTables.fixedHeader.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js"></script>
</head>
<body>
    <div class="container">
        <h2>{report_title}</h2>
        
        <div class="info-box">
            <div style="text-align: center; margin-bottom: 15px; min-height: 80px;">
                {category_svg}
            </div>
            <p><strong>{category_intro_label}</strong> {category_description}</p>
            <p>
                For more details on categories, visit here: <a href="https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-isoform-classification:-categories-and-subcategories" target="_blank">SQANTI3 Wiki: Categories</a>.
                <br>
                For a detailed explanation of each column of the table (SQANTI3 classification), see here: <a href="https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC#classifcols" target="_blank">Column Glossary</a>.
            </p>
            <hr style="border-top: 1px solid #dee2e6; margin: 15px 0;">
            <p><strong>How to Filter:</strong></p>
            <ul style="margin-bottom: 5px; padding-left: 20px;">
                <li><strong>Numeric Columns (e.g., length, exons):</strong> Use ranges like <code>100:1000</code> (values between 100 and 1000), <code>100:</code> (>100), or <code>:1000</code> (<1000). Single numbers perform an exact match.</li>
                <li><strong>Dropdowns:</strong> Select a value to filter for an exact match.</li>
                <li><strong>Text Columns:</strong> Type to search for substrings (case-insensitive).</li>
            </ul>
            <hr style="border-top: 1px solid #dee2e6; margin: 15px 0;">
            <p><strong>Generate Trix String (for UCSC Genome Browser search):</strong></p>
            <ul style="margin-bottom: 5px; padding-left: 20px;">
                <li>Use dropdown filters and click "Generate Trix String" to get a search based on your filter criteria.</li>
                <li><em>Note:</em> Range filters (e.g., <code>100:1000</code>) are not supported by Trix and will be ignored, only exact numeric matches are supported.</li>
                <li>Search terms automatically use underscores (e.g., <code>structural_category_full_splice_match</code>, <code>strand_plus</code>, <code>coding_coding</code>).</li>
                <li>Click on a row to select it (highlighted in blue), then click "Generate Trix String" to get search terms for that specific isoform.</li>
            </ul>
        </div>

        <p>Total transcripts: {count} | <a href="javascript:window.close();">Close Window</a></p>
        
        <table id="isoformTable" class="display stripe hover row-border order-column" style="width:100%">
            <thead>
                <tr>
                    {headers}
                </tr>
                <tr class="filters">
                    {headers}
                </tr>
            </thead>
        </table>
    </div>

    <script>
        // The data is embedded directly into the HTML to make it self-contained (single file)
        var tableData = {data_json};
        var tableColumns = {columns_json};
        var numericColumns = {numeric_columns_json};
        var categoricalColumns = {categorical_columns_json};
        var currentCategory = "{category_safe}";
        var includeCategoryFilter = {include_category_filter};

        function normalizeForTrix(value) {{
            if (value === undefined || value === null) {{
                return '';
            }}
            var normalized = String(value).trim().toLowerCase();
            if (normalized === '+') {{
                return 'plus';
            }}
            if (normalized === '-') {{
                return 'minus';
            }}
            normalized = normalized.replace(/[^\w]+/g, '_');
            normalized = normalized.replace(/_+/g, '_');
            normalized = normalized.replace(/^_+|_+$/g, '');
            return normalized;
        }}

        // Custom filtering function for range search and special exact match cases
        $.fn.dataTable.ext.search.push(
            function( settings, data, dataIndex ) {{
                var result = true;
                
                // Iterate over all columns to check for active range filters
                // We use settings.aoColumns to get column info
                for (var i = 0; i < settings.aoColumns.length; i++) {{
                    var colName = settings.aoColumns[i].data;
                    
                    // Find the input value for this column
                    // We target the input in the visible scroll header first
                    var input = $('.dataTables_scrollHead .filters [data-index="'+i+'"]');
                    if (input.length === 0) {{
                        // Fallback to any input with this index (e.g. if not scrolling)
                        input = $('.filters [data-index="'+i+'"]');
                    }}
                    
                    var filterVal = input.val();
                    if (!filterVal) continue;
                    
                    // Special case: Exact match for categorical columns (dropdowns)
                    if (categoricalColumns.includes(colName)) {{
                        // Exact string match (case sensitive)
                        // We convert data to string just in case (e.g. numbers treated as categories)
                        if (String(data[i]) !== filterVal) {{
                            result = false;
                            break;
                        }}
                        continue; 
                    }}
                    
                    // Only check numeric columns for range logic
                    if (!numericColumns.includes(colName)) continue;
                    
                    // Only apply custom range logic if filter contains ':'
                    if (filterVal.indexOf(':') !== -1) {{
                        var parts = filterVal.split(':');
                        var minStr = parts[0].trim();
                        var maxStr = parts[1].trim();
                        
                        var min = minStr !== "" ? parseFloat(minStr) : -Infinity;
                        var max = maxStr !== "" ? parseFloat(maxStr) : Infinity;
                        
                        var cellVal = parseFloat(data[i]);
                        
                        if (isNaN(cellVal)) {{
                            result = false;
                            break;
                        }}
                        
                        if (cellVal < min || cellVal > max) {{
                            result = false;
                            break;
                        }}
                    }} else {{
                        // Exact numeric match (not substring match)
                        var searchNum = parseFloat(filterVal);
                        var cellVal = parseFloat(data[i]);
                        
                        // If the filter value is a valid number
                        if (!isNaN(searchNum)) {{
                            // If cell is not a number (e.g. empty), it doesn't match
                            if (isNaN(cellVal)) {{
                                result = false;
                                break;
                            }}
                            
                            // Exact match check
                            if (cellVal !== searchNum) {{
                                result = false;
                                break;
                            }}
                        }}
                        // If filter is not a number (e.g. text search in numeric col)
                        else if (filterVal === data[i]) {{
                             // Exact string match allowed (e.g. for 'NA')
                        }}
                        else {{
                            // Filter is text but doesn't match exact string
                             result = false; 
                             break;
                        }}
                    }}
                }}
                
                return result;
            }}
        );

        $(document).ready(function() {{
            // Setup - add inputs or dropdowns
            $('#isoformTable thead tr.filters th').each( function (i) {{
                var title = $(this).text();
                
                if (categoricalColumns.includes(title)) {{
                    // Create dropdown
                    var select = $('<select data-index="'+i+'"><option value="">All</option></select>');
                    
                    // Get unique values from tableData for this column
                    // We use Set to get unique values
                    var unique = [...new Set(tableData.map(item => item[title]))].sort();
                    
                    unique.forEach(function(d) {{
                        // Skip empty values for options
                        if (d !== "" && d !== null && d !== undefined) {{
                            select.append('<option value="'+d+'">'+d+'</option>');
                        }}
                    }});
                    
                    $(this).html(select);
                }} else {{
                    // Create text input
                    $(this).html( '<input type="text" placeholder="Search '+title+'" data-index="'+i+'" />' );
                }}
            }} );
         
            // Initialize DataTable
            var table = $('#isoformTable').DataTable({{
                data: tableData,
                columns: tableColumns,
                dom: 'Bfrtip',
                pageLength: 25,
                lengthMenu: [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]],
                scrollX: true,     // Enable horizontal scrolling
                scrollY: '70vh',   // Enable vertical scrolling with fixed header
                scrollCollapse: true,
                autoWidth: true,
                orderCellsTop: true, // Important: tells DataTables the first row is for ordering
                buttons: [
                    'pageLength',
                    {{
                        text: 'Generate Trix String',
                        action: function ( e, dt, node, config ) {{
                            // Check if a row is selected
                            var selectedData = dt.rows('.selected').data();
                            
                            if (selectedData.length > 0) {{
                                // Generate Trix string from selected row
                                var rowData = selectedData[0];
                                var queryParts = [];
                                
                                // Include structural category context if available
                                if (includeCategoryFilter) {{
                                    var normalizedCategory = normalizeForTrix(currentCategory);
                                    if (normalizedCategory) {{
                                        queryParts.push('structural_category_' + normalizedCategory);
                                    }}
                                }} else if (rowData['structural_category']) {{
                                    var normalizedCategory = normalizeForTrix(rowData['structural_category']);
                                    if (normalizedCategory) {{
                                        queryParts.push('structural_category_' + normalizedCategory);
                                    }}
                                }}
                                
                                // Columns to include in the Trix string (categorical/useful for searching)
                                // Note: structural_category is already added above
                                var trixCols = ['subcategory', 'strand', 'coding', 
                                                'associated_gene', 'FSM_class', 'RTS_stage', 'all_canonical',
                                                'bite', 'predicted_NMD', 'within_CAGE_peak', 'polyA_motif_found'];
                                
                                trixCols.forEach(function(col) {{
                                    var value = rowData[col];
                                    if (value !== undefined && value !== '' && value !== null) {{
                                        var normalizedCol = normalizeForTrix(col);
                                        var normalizedVal = normalizeForTrix(value);
                                        if (normalizedCol && normalizedVal) {{
                                            queryParts.push(normalizedCol + '_' + normalizedVal);
                                        }}
                                    }}
                                }});
                                
                                // Also add the isoform name for direct search
                                var isoformCol = 'isoform';
                                if (rowData[isoformCol]) {{
                                    queryParts.unshift(rowData[isoformCol]);
                                }}
                                
                                var uniqueParts = [];
                                var seenParts = new Set();
                                queryParts.forEach(function(part) {{
                                    if (part && !seenParts.has(part)) {{
                                        seenParts.add(part);
                                        uniqueParts.push(part);
                                    }}
                                }});

                                var searchString = uniqueParts.join(" ");
                                prompt("Trix search string for selected isoform:\\n(You can use just the isoform ID, or any combination of the terms below)", searchString);
                            }} else {{
                                // Generate from active filters (non-range values only)
                                var queryParts = [];
                                var api = dt;
                                var settings = api.settings()[0];
                                
                                // Include structural category context if available
                                if (includeCategoryFilter) {{
                                    var normalizedCategory = normalizeForTrix(currentCategory);
                                    if (normalizedCategory) {{
                                        queryParts.push('structural_category_' + normalizedCategory);
                                    }}
                                }}
                                
                                api.columns().every(function(i) {{
                                    var colName = settings.aoColumns[i].data;
                                    // Skip structural_category since we already added it
                                    if (includeCategoryFilter && colName === 'structural_category') return;
                                    
                                    var input = $('.dataTables_scrollHead .filters [data-index="'+i+'"]');
                                    if (input.length === 0) {{
                                        input = $('.filters [data-index="'+i+'"]');
                                    }}
                                    
                                    var val = input.val();
                                    if (val) {{
                                        // Skip range filters (contain :) - Trix doesn't support ranges
                                        if (val.indexOf(':') !== -1) {{
                                            return; // Skip this column
                                        }}
                                        var normalizedCol = normalizeForTrix(colName);
                                        var normalizedVal = normalizeForTrix(val);
                                        if (normalizedCol && normalizedVal) {{
                                            queryParts.push(normalizedCol + '_' + normalizedVal);
                                        }}
                                    }}
                                }});
                                
                                var uniqueParts = [];
                                var seenParts = new Set();
                                queryParts.forEach(function(part) {{
                                    if (part && !seenParts.has(part)) {{
                                        seenParts.add(part);
                                        uniqueParts.push(part);
                                    }}
                                }});

                                // queryParts always has at least structural_category now
                                prompt("Trix search string from filters:\\n(Note: Range filters are ignored - Trix doesn't support ranges)", uniqueParts.join(" "));
                            }}
                        }}
                    }},
                    {{
                        extend: 'excelHtml5',
                        title: '{export_title}'
                    }},
                    {{
                        extend: 'csvHtml5',
                        title: '{export_title}'
                    }},
                    'copy'
                ],
                initComplete: function () {{
                    // No manual event binding here anymore, we use delegated events below
                }}
            }});
            
            // Row click handler for selection (toggle)
            $('#isoformTable tbody').on('click', 'tr', function () {{
                if ($(this).hasClass('selected')) {{
                    $(this).removeClass('selected');
                }} else {{
                    table.$('tr.selected').removeClass('selected');
                    $(this).addClass('selected');
                }}
            }});
            
            // Delegated event listener for inputs AND selects
            $(table.table().container()).on('keyup change clear', 'thead input, thead select', function () {{
                var input = $(this);
                var colIdx = input.data('index');
                var val = this.value;
                
                // Sync value to all inputs/selects for this column
                var tagName = input.prop("tagName");
                $(tagName + '[data-index="'+colIdx+'"]').val(val);
                
                var column = table.column(colIdx);
                var colName = column.dataSrc();
                
                // If numeric column OR categorical column, ALWAYS use custom search (clear standard search)
                if (numericColumns.includes(colName) || categoricalColumns.includes(colName)) {{
                    if (column.search() !== "") {{
                        column.search(""); 
                    }}
                    table.draw(); // Trigger custom filter
                }} else {{
                    // Standard search behavior for other columns
                    if (column.search() !== val) {{
                        column.search(val).draw();
                    }}
                }}
            }});
        }});
    </script>
</body>
</html>
"""

# Definitions from SQANTI3 Wiki
CATEGORY_DEFINITIONS = {
    "full-splice_match": "<b>FSM (Full Splice Match):</b> The reference and query isoform have the same number of exons and each internal junction agree. The exact 5' start and 3' end can differ by any amount.",
    "incomplete-splice_match": "<b>ISM (Incomplete Splice Match):</b> The query isoform has fewer 5' exons than the reference, but each internal junction agree. The exact 5' start and 3' end can differ by any amount.",
    "novel_in_catalog": "<b>NIC (Novel In Catalog):</b> The query isoform does not have a FSM or ISM match, but is using a combination of known donor/acceptor sites.",
    "novel_not_in_catalog": "<b>NNC (Novel Not in Catalog):</b> The query isoform does not have a FSM or ISM match, and has at least one donor or acceptor site that is not annotated.",
    "antisense": "<b>Antisense:</b> The query isoform does not overlap a same-strand reference gene but is anti-sense to an annotated gene.",
    "genic_intron": "<b>Genic Intron:</b> The query isoform is completely contained within an annotated intron.",
    "genic": "<b>Genic Genomic:</b> The query isoform overlaps with introns and exons.",
    "intergenic": "<b>Intergenic:</b> The query isoform is in the intergenic region.",
    "fusion": "<b>Fusion:</b> The query isoform spans two or more annotated genes.",
    "NA": "No structural category assigned."
}

# Preferred display order for structural categories
CATEGORY_ORDER = [
    "full-splice_match",
    "incomplete-splice_match",
    "novel_in_catalog",
    "novel_not_in_catalog",
    "genic",
    "antisense",
    "fusion",
    "intergenic",
    "genic_intron",
    "NA"
]

# SVG Diagrams for structural categories - Text labels moved to the right side
CATEGORY_SVGS = {
    "full-splice_match": '''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="190" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="80" y="8" width="40" height="8" fill="black" rx="2"/><rect x="150" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><line x1="10" y1="52" x2="190" y2="52" stroke="#6BAED6" stroke-width="2" /><rect x="10" y="48" width="40" height="8" fill="#6BAED6" rx="2"/><rect x="80" y="48" width="40" height="8" fill="#6BAED6" rx="2"/><rect x="150" y="48" width="40" height="8" fill="#6BAED6" rx="2"/><text x="220" y="55" font-family="sans-serif" font-size="12" fill="#6BAED6" font-weight="bold">Isoform (FSM)</text></svg>''',
    
    "incomplete-splice_match": '''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="190" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="80" y="8" width="40" height="8" fill="black" rx="2"/><rect x="150" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><line x1="80" y1="52" x2="190" y2="52" stroke="#FC8D59" stroke-width="2" /><rect x="80" y="48" width="40" height="8" fill="#FC8D59" rx="2"/><rect x="150" y="48" width="40" height="8" fill="#FC8D59" rx="2"/><text x="220" y="55" font-family="sans-serif" font-size="12" fill="#FC8D59" font-weight="bold">Isoform (ISM)</text></svg>''',
    
    "novel_in_catalog": '''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="190" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="80" y="8" width="40" height="8" fill="black" rx="2"/><rect x="150" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><line x1="10" y1="52" x2="190" y2="52" stroke="#78C679" stroke-width="2" /><rect x="10" y="48" width="40" height="8" fill="#78C679" rx="2"/><rect x="150" y="48" width="40" height="8" fill="#78C679" rx="2"/><text x="220" y="55" font-family="sans-serif" font-size="12" fill="#78C679" font-weight="bold">Isoform (NIC)</text></svg>''',
    
    "novel_not_in_catalog": '''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="190" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="80" y="8" width="40" height="8" fill="black" rx="2"/><rect x="150" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><line x1="10" y1="52" x2="190" y2="52" stroke="#D62F4B" stroke-width="2" /><rect x="10" y="48" width="40" height="8" fill="#D62F4B" rx="2"/><rect x="80" y="48" width="60" height="8" fill="#D62F4B" rx="2"/><rect x="150" y="48" width="40" height="8" fill="#D62F4B" rx="2"/><text x="130" y="40" font-size="10" fill="#D62F4B" text-anchor="middle">New Site</text><text x="220" y="55" font-family="sans-serif" font-size="12" fill="#D62F4B" font-weight="bold">Isoform (NNC)</text></svg>''',
    
    "genic_intron": '''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="190" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="150" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><rect x="80" y="48" width="40" height="8" fill="#41B6C4" rx="2"/><text x="220" y="55" font-family="sans-serif" font-size="12" fill="#41B6C4" font-weight="bold">Isoform (Intron)</text></svg>''',
    
    "genic": '''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="120" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="80" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><rect x="40" y="48" width="60" height="8" fill="#969696" rx="2"/><text x="70" y="70" font-size="10" fill="#969696" text-anchor="middle">Overlap</text><text x="220" y="55" font-family="sans-serif" font-size="12" fill="#969696" font-weight="bold">Isoform (Genic)</text></svg>''',
    
    "antisense": '''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="12" x2="120" y2="12" stroke="black" stroke-width="2" /><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><rect x="80" y="8" width="40" height="8" fill="black" rx="2"/><text x="130" y="15" font-size="12">→</text><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference (+)</text><line x1="10" y1="52" x2="120" y2="52" stroke="#66C2A4" stroke-width="2" /><rect x="10" y="48" width="40" height="8" fill="#66C2A4" rx="2"/><rect x="80" y="48" width="40" height="8" fill="#66C2A4" rx="2"/><text x="130" y="55" font-size="12" fill="#66C2A4">←</text><text x="220" y="55" font-family="sans-serif" font-size="12" fill="#66C2A4" font-weight="bold">Isoform (-)</text></svg>''',
    
    "intergenic": '''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><rect x="140" y="48" width="40" height="8" fill="#E9967A" rx="2"/><text x="220" y="55" font-family="sans-serif" font-size="12" fill="#E9967A" font-weight="bold">Isoform (Inter)</text></svg>''',
    
    "fusion": '''<svg width="550" height="110" viewBox="0 0 400 80" xmlns="http://www.w3.org/2000/svg"><rect x="10" y="8" width="40" height="8" fill="black" rx="2"/><text x="55" y="15" font-size="10">Gene A</text><rect x="140" y="8" width="40" height="8" fill="black" rx="2"/><text x="185" y="15" font-size="10">Gene B</text><text x="220" y="15" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text><line x1="10" y1="52" x2="180" y2="52" stroke="#DAA520" stroke-width="2" /><rect x="10" y="48" width="40" height="8" fill="#DAA520" rx="2"/><rect x="140" y="48" width="40" height="8" fill="#DAA520" rx="2"/><text x="220" y="55" font-family="sans-serif" font-size="12" fill="#DAA520" font-weight="bold">Isoform (Fusion)</text></svg>'''
}

REFERENCE_SVG = '''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="30" x2="190" y2="30" stroke="black" stroke-width="2" /><rect x="10" y="26" width="40" height="8" fill="black" rx="2"/><rect x="80" y="26" width="40" height="8" fill="black" rx="2"/><rect x="150" y="26" width="40" height="8" fill="black" rx="2"/><text x="220" y="33" font-family="sans-serif" font-size="12" fill="black" font-weight="bold">Reference</text></svg>'''

# Isoform-only SVGs for complete transcriptome overview (no repeated reference)
CATEGORY_SVGS_OVERVIEW = {
    "full-splice_match": '''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="30" x2="190" y2="30" stroke="#6BAED6" stroke-width="2" /><rect x="10" y="26" width="40" height="8" fill="#6BAED6" rx="2"/><rect x="80" y="26" width="40" height="8" fill="#6BAED6" rx="2"/><rect x="150" y="26" width="40" height="8" fill="#6BAED6" rx="2"/><text x="220" y="33" font-family="sans-serif" font-size="12" fill="#6BAED6" font-weight="bold">Isoform (FSM)</text></svg>''',
    "incomplete-splice_match": '''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><line x1="80" y1="30" x2="190" y2="30" stroke="#FC8D59" stroke-width="2" /><rect x="80" y="26" width="40" height="8" fill="#FC8D59" rx="2"/><rect x="150" y="26" width="40" height="8" fill="#FC8D59" rx="2"/><text x="220" y="33" font-family="sans-serif" font-size="12" fill="#FC8D59" font-weight="bold">Isoform (ISM)</text></svg>''',
    "novel_in_catalog": '''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="30" x2="190" y2="30" stroke="#78C679" stroke-width="2" /><rect x="10" y="26" width="40" height="8" fill="#78C679" rx="2"/><rect x="150" y="26" width="40" height="8" fill="#78C679" rx="2"/><text x="220" y="33" font-family="sans-serif" font-size="12" fill="#78C679" font-weight="bold">Isoform (NIC)</text></svg>''',
    "novel_not_in_catalog": '''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="30" x2="190" y2="30" stroke="#D62F4B" stroke-width="2" /><rect x="10" y="26" width="40" height="8" fill="#D62F4B" rx="2"/><rect x="80" y="26" width="60" height="8" fill="#D62F4B" rx="2"/><rect x="150" y="26" width="40" height="8" fill="#D62F4B" rx="2"/><text x="130" y="18" font-size="10" fill="#D62F4B" text-anchor="middle">New Site</text><text x="220" y="33" font-family="sans-serif" font-size="12" fill="#D62F4B" font-weight="bold">Isoform (NNC)</text></svg>''',
    "genic_intron": '''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><rect x="80" y="26" width="40" height="8" fill="#41B6C4" rx="2"/><text x="220" y="33" font-family="sans-serif" font-size="12" fill="#41B6C4" font-weight="bold">Isoform (Intron)</text></svg>''',
    "genic": '''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><rect x="40" y="26" width="60" height="8" fill="#969696" rx="2"/><text x="70" y="48" font-size="10" fill="#969696" text-anchor="middle">Overlap</text><text x="220" y="33" font-family="sans-serif" font-size="12" fill="#969696" font-weight="bold">Isoform (Genic)</text></svg>''',
    "antisense": '''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="30" x2="120" y2="30" stroke="#66C2A4" stroke-width="2" /><rect x="10" y="26" width="40" height="8" fill="#66C2A4" rx="2"/><rect x="80" y="26" width="40" height="8" fill="#66C2A4" rx="2"/><text x="130" y="33" font-size="12" fill="#66C2A4">←</text><text x="220" y="33" font-family="sans-serif" font-size="12" fill="#66C2A4" font-weight="bold">Isoform (-)</text></svg>''',
    "intergenic": '''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><rect x="140" y="26" width="40" height="8" fill="#E9967A" rx="2"/><text x="220" y="33" font-family="sans-serif" font-size="12" fill="#E9967A" font-weight="bold">Isoform (Inter)</text></svg>''',
    "fusion": '''<svg width="550" height="70" viewBox="0 10 400 50" xmlns="http://www.w3.org/2000/svg"><line x1="10" y1="30" x2="180" y2="30" stroke="#DAA520" stroke-width="2" /><rect x="10" y="26" width="40" height="8" fill="#DAA520" rx="2"/><rect x="140" y="26" width="40" height="8" fill="#DAA520" rx="2"/><text x="220" y="33" font-family="sans-serif" font-size="12" fill="#DAA520" font-weight="bold">Isoform (Fusion)</text></svg>'''
}

def build_category_overview_html(categories_present=None, palette=None):
    """
    Build an HTML block with SVGs and definitions for structural categories.
    
    Args:
        categories_present: Set of categories to include (defaults to all)
        palette: Optional custom color palette (dict of category -> RGB tuple)
    """
    if categories_present is None:
        categories_present = set(CATEGORY_ORDER)
    else:
        categories_present = set(categories_present)
    
    blocks = [
        f"""
        <div style="display: flex; align-items: center; gap: 8px; margin: 2px 0; padding: 0;">
            <div style="flex: 0 0 260px; text-align: center; line-height: 0;">{REFERENCE_SVG}</div>
            <div style="flex: 1;"><strong>Reference:</strong> Shared reference drawing used for all categories.</div>
        </div>
        """.strip()
    ]
    for cat in CATEGORY_ORDER:
        if cat not in categories_present:
            continue
        definition = CATEGORY_DEFINITIONS.get(cat)
        if not definition:
            continue
        # Generate SVG with custom colors if palette provided
        svg = generate_category_svg(cat, palette, include_reference=False)
        blocks.append(
            f"""
            <div style="display: flex; align-items: center; gap: 8px; margin: 2px 0; padding: 0;">
                <div style="flex: 0 0 260px; text-align: center; line-height: 0;">{svg}</div>
                <div style="flex: 1;">{definition}</div>
            </div>
            """.strip()
        )
    
    if not blocks:
        return "<p>No structural category definitions available.</p>"
    
    return "\n".join(blocks)

def generate_html_reports(classification_file, output_dir, include_sequences=False, custom_palette=None):
    """
    Generate HTML reports for each structural category in the classification file.
    
    Args:
        classification_file: Path to SQANTI3 classification file
        output_dir: Output directory for HTML reports
        include_sequences: Whether to include ORF_seq column
        custom_palette: Optional custom color palette (dict of category -> RGB tuple)
    """
    # Opt-in to future pandas behavior to silence downcasting warnings
    try:
        pd.set_option('future.no_silent_downcasting', True)
    except Exception:
        pass

    logger.info(f"Reading classification file: {classification_file}")
    try:
        df = pd.read_csv(classification_file, sep='\t')
    except Exception as e:
        logger.exception(f"Error reading classification file: {e}")
        sys.exit(1)

    # Clean column names (remove spaces/special chars for JSON compatibility if needed, 
    # but DataTables handles strings fine. We'll keep original names for display)
    
    # Optionally exclude ORF_seq column
    if not include_sequences and 'ORF_seq' in df.columns:
        logger.info("Excluding 'ORF_seq' column (use --include-sequences to keep it)")
        df = df.drop(columns=['ORF_seq'])

    if 'structural_category' not in df.columns:
        logger.error("Error: 'structural_category' column not found in classification file.")
        sys.exit(1)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    categories = df['structural_category'].unique()
    logger.info(f"Found {len(categories)} structural categories.")

    # Columns that should have dropdown filters
    categorical_cols = [
        'strand', 'structural_category', 'subcategory', 'RTS_stage', 
        'all_canonical', 'min_cov_pos', 'bite', 'FSM_class', 'coding', 
        'predicted_NMD', 'within_CAGE_peak', 'polyA_motif_found'
    ]

    # Full transcriptome report (all categories)
    full_df = df.copy()
    numeric_cols_all = full_df.select_dtypes(include=['number']).columns.tolist()
    full_df = full_df.fillna('')
    
    columns_json = []
    headers_html = ""
    for col in full_df.columns:
        columns_json.append({"data": col, "title": col})
        headers_html += f"<th>{col}</th>"
    
    data_json = full_df.to_json(orient='records')
    category_overview = build_category_overview_html(categories_present=categories, palette=custom_palette)
    full_report_title = "SQANTI3 Analysis: complete transcriptome"
    full_export_title = "SQANTI3_complete_transcriptome_Isoforms"
    
    full_filename = "complete_transcriptome_isoforms.html"
    full_category_safe = normalize_trix_token("complete transcriptome") or "complete_transcriptome"
    full_path = output_path / full_filename
    full_html = HTML_TEMPLATE.format(
        category="complete transcriptome",
        category_safe=full_category_safe,
        report_title=full_report_title,
        export_title=full_export_title,
        category_intro_label="Structural Categories:",
        category_description="All structural categories and definitions are shown below.",
        category_svg=category_overview,
        count=len(full_df),
        headers=headers_html,
        data_json=data_json,
        columns_json=json.dumps(columns_json),
        numeric_columns_json=json.dumps(numeric_cols_all),
        categorical_columns_json=json.dumps(categorical_cols),
        include_category_filter="false"
    )
    
    with open(full_path, 'w') as f:
        f.write(full_html)
    logger.info(f"Generated full transcriptome report: {full_path}")

    for category in categories:
        if pd.isna(category) or category == 'NA':
            continue

        cat_str = str(category)
        logger.info(f"Processing category: {cat_str}")
        cat_safe = normalize_trix_token(cat_str) or cat_str.replace(' ', '_')
        
        # Get category definition
        cat_def = CATEGORY_DEFINITIONS.get(cat_str, CATEGORY_DEFINITIONS.get(cat_str.replace(' ', '_'), "No definition available."))
        
        # Get category SVG (with custom colors if palette provided)
        cat_svg = generate_category_svg(cat_str, custom_palette, include_reference=True)
        if not cat_svg:
            cat_svg = generate_category_svg(cat_str.replace(' ', '_'), custom_palette, include_reference=True)
        
        # Filter data
        cat_df = df[df['structural_category'] == category].copy()
        
        # Identify numeric columns for range filtering (before converting to string/filling NaNs)
        numeric_cols = cat_df.select_dtypes(include=['number']).columns.tolist()
        
        # Replace NaNs with empty string for JSON/HTML display
        cat_df = cat_df.fillna('')
        
        # Prepare data for DataTables
        # columns definition: [{ data: 'col_name', title: 'Col Name' }, ...]
        columns_json = []
        headers_html = ""
        
        for col in cat_df.columns:
            columns_json.append({"data": col, "title": col})
            headers_html += f"<th>{col}</th>"

        # Convert data to JSON
        data_json = cat_df.to_json(orient='records')

        # Safe filename
        safe_cat = cat_str.replace(' ', '_').replace('/', '_')
        filename = f"{safe_cat}_isoforms.html"
        file_path = output_path / filename

        # Fill template
        report_title = f"SQANTI3 Analysis: {cat_str} Isoforms"
        export_title = f"SQANTI3_{safe_cat}_Isoforms"
        html_content = HTML_TEMPLATE.format(
            category=cat_str,
            category_safe=cat_safe,  # Used for Trix string generation (matches Trix index format)
            report_title=report_title,
            export_title=export_title,
            category_intro_label="Category Definition:",
            category_description=cat_def,
            category_svg=cat_svg,
            count=len(cat_df),
            headers=headers_html,
            data_json=data_json,
            columns_json=json.dumps(columns_json),
            numeric_columns_json=json.dumps(numeric_cols),
            categorical_columns_json=json.dumps(categorical_cols),
            include_category_filter="true"
        )

        with open(file_path, 'w') as f:
            f.write(html_content)
        
        logger.info(f"  -> Generated {file_path}")

    logger.info("Done! All reports generated.")

def main():
    parser = argparse.ArgumentParser(description="Generate filterable HTML tables for SQANTI3 isoform categories.")
    parser.add_argument("--classification", required=True, help="Path to SQANTI3 classification file")
    parser.add_argument("--output-dir", default="isoform_reports", help="Directory to save HTML files")
    parser.add_argument("--include-sequences", action="store_true", help="Include ORF_seq column (warning: large file size)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.classification):
        logger.error(f"Error: File not found: {args.classification}")
        sys.exit(1)

    generate_html_reports(args.classification, args.output_dir, args.include_sequences)

if __name__ == "__main__":
    main()
