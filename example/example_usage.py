#!/usr/bin/env python3
"""
Example usage script for SQANTI3 to UCSC Genome Browser integration

This script demonstrates how to use the integration tool with example data.

Author: Carolina Monzo
"""

import os
import sys
import subprocess
from pathlib import Path

def run_example():
    """Run an example conversion"""
    project_root = Path(__file__).parent.parent

    print("SQANTI3 to UCSC Genome Browser Integration - Example Usage")
    print("=" * 60)

    # Check if the main module exists
    main_module = "sqanti_browser"
    main_script = "sqanti_browser.py"  # for existence check
    filter_script = "src/filter_isoforms.py"

    if not (project_root / main_script).exists():
        print(f"Error: {main_script} not found in project root")
        return False

    if not (project_root / filter_script).exists():
        print(f"Warning: {filter_script} not found in project")

    # Example parameters (paths relative to project root)
    example_gtf = "example/SQANTI3_QC_output/example_corrected.gtf"
    example_classification = "example/SQANTI3_QC_output/example_classification.txt"
    output_dir = "example_output"
    genome = "hg38"
    
    print(f"Example parameters:")
    print(f"  GTF file: {example_gtf}")
    print(f"  Classification file: {example_classification}")
    print(f"  Output directory: {output_dir}")
    print(f"  Genome: {genome}")
    print()
    
    # Check if example files exist
    if not (project_root / example_gtf).exists():
        print(f"Note: Example GTF file not found: {example_gtf}")
        print("This is expected if you haven't downloaded the SQANTI3 example data yet.")
        print()
        print("To download example data:")
        print("  git clone https://github.com/ConesaLab/SQANTI3.git")
        print("  cp -r SQANTI3/example/SQANTI3_QC_output ./example/")
        print()
    
    if not (project_root / example_classification).exists():
        print(f"Note: Example classification file not found: {example_classification}")
        print("This is expected if you haven't downloaded the SQANTI3 example data yet.")
        print()
    
    # Show the command that would be run
    print("To run the conversion with your own data, use:")
    print()
    print(f"python -m {main_module} \\")
    print(f"    --gtf <your_corrected.gtf> \\")
    print(f"    --classification <your_classification.txt> \\")
    print(f"    --output <output_directory> \\")
    print(f"    --genome <genome_assembly>")
    print()
    print("For GitHub integration (recommended for UCSC Genome Browser):")
    print(f"python -m {main_module} \\")
    print(f"    --gtf <your_corrected.gtf> \\")
    print(f"    --classification <your_classification.txt> \\")
    print(f"    --output <output_directory> \\")
    print(f"    --genome <genome_assembly>")
    print()
    print("Optional enhancements:")
    print("- Generate interactive HTML tables: --tables")
    print("- Add STAR splice junctions track from SJ.out.tab: --star-sj <path>")
    print("- Add CAGE peaks for TSS validation: --CAGE-peak <CAGE_peaks.bed>")
    print("- Add PolyA peaks for TTS validation: --polyA-peak <polyA_peaks.bed>")
    print("- Add reference annotation for comparison: --refGTF <reference.gtf>")
    print("  (Validation tracks show by default when the hub is opened; colors are embedded in bigBed)")
    print("- Use a genome .2bit to compute chrom.sizes: --twobit <genome.2bit>")
    print("- Sort isoforms: --sort-by <metric> (default none = preserve original GTF order). Options: FL, iso_exp, length, diff_to_TSS, etc.")
    print("- Skip category-specific tracks: --no-category-tracks")
    print("- Only specific category tracks: --category-tracks FSM,ISM,NIC (abbrev: FSM, ISM, NIC, NNC, antisense, genic_intron, genic_genomic, intergenic, fusion)")
    print("- Disable highlighted top FL isoforms: --no-highlight")
    print("- Custom color palette: --my-palette <palette.json> (see example/example_palette.json)")
    print("- Validate tools/inputs only (no outputs): --validate-only")
    print("- Dry run (build enhanced BED, skip bigBed/hub): --dry-run")
    print()
    print("Note: For info on obtaining CAGE and PolyA peak files, see the SQANTI3 wiki:")
    print("      https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control")
    print()
    print("Note: For the example reference GTF (chr22), download from:")
    print("      https://github.com/ConesaLab/SQANTI3/blob/master/data/reference/gencode.v38.basic_chr22.gtf")
    print()
    print("With --my-palette, track HTML pages and table reports use your custom colors.")
    print("Trix search uses underscore format (e.g., structural_category_full_splice_match).")
    print("To view specific isoforms: use HTML tables (--tables), select rows, Generate Filter String,")
    print("  then UCSC Tools -> Table Browser -> Identifiers -> paste list -> Output format: custom track.")
    print()
    print("Generate interactive HTML reports for exploration:")
    print(f"python {filter_script} \\")
    print(f"    --classification <your_classification.txt> \\")
    print(f"    --output-dir <report_output_dir>")
    print()
    print("Run with SQANTI3 example dataset (after copying example data):")
    print(f"python -m {main_module} \\")
    print(f"    --gtf {example_gtf} \\")
    print(f"    --classification {example_classification} \\")
    print(f"    --output {output_dir} \\")
    print(f"    --genome {genome}")
    print()
    
    # Show help
    print("For detailed help and all options:")
    print(f"python -m {main_module} --help")
    print()
    
    # Check if we can run the script (python -m ensures active env is used)
    try:
        result = subprocess.run([sys.executable, "-m", main_module, "--help"], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            print("✓ Main script is working correctly")
        else:
            print("✗ Main script help command failed")
            return False
    except Exception as e:
        print(f"✗ Could not run main script: {e}")
        return False
    
    print()
    print("Example workflow:")
    print("1. Run SQANTI3 on your transcriptome data")
    print("2. Use this tool to convert the output to bigBed format")
    print("3. (Optional) Generate interactive HTML reports for detailed exploration")
    print("4. Upload the generated hub files to a web server")
    print("5. Add the hub to UCSC Genome Browser")
    print("6. View and filter your transcripts by structural category")
    print()
    print("For detailed instructions, see README.md")
    
    return True

def show_file_structure():
    """Show the expected file structure"""
    print("Expected file structure after conversion:")
    print()
    print("output_directory/")
    print("├── hub.txt                              # Hub configuration")
    print("├── genomes.txt                          # Genome mapping")
    print("├── README.md                            # Hub documentation")
    print("├── hg38/                                # Genome-specific directory")
    print("│   ├── hg38_sqanti3.bb                  # Main bigBed (all transcripts)")
    print("│   ├── hg38_sqanti3_full-splice_match.bb    # Category-specific bigBed")
    print("│   ├── hg38_sqanti3_novel_in_catalog.bb     # Category-specific bigBed")
    print("│   ├── hg38_sqanti3_*.bb                    # Other category bigBeds")
    print("│   ├── hg38_star_sj.bb                  # STAR junctions (optional)")
    print("│   ├── hg38_cage_peaks.bb               # CAGE peaks (optional)")
    print("│   ├── hg38_polya_peaks.bb              # PolyA peaks (optional)")
    print("│   ├── hg38_reference.bb                # Reference annotation (optional)")
    print("│   ├── hg38_sqanti3_track.html              # Main track documentation")
    print("│   ├── hg38_sqanti3_*.html                  # Category track documentation")
    print("│   ├── hg38_star_sj_track.html              # STAR track doc (optional)")
    print("│   ├── hg38_cage_peaks_track.html           # CAGE track doc (optional)")
    print("│   ├── hg38_polya_peaks_track.html          # PolyA track doc (optional)")
    print("│   ├── hg38_reference_track.html            # Reference track doc (optional)")
    print("│   ├── trackDb.txt                      # Track configuration")
    print("│   ├── groups.txt                       # Track groups")
    print("│   ├── trix.ix                          # Trix search index")
    print("│   └── trix.ixx                         # Trix search index")
    print("└── table_reports/                       # (if --tables flag used)")
    print("    ├── complete_transcriptome_isoforms.html  # All categories combined")
    print("    ├── full-splice_match_isoforms.html  # Interactive HTML table")
    print("    └── ...                              # One HTML per category")
    print()

def main():
    """Main function"""
    if not run_example():
        sys.exit(1)
    
    print()
    show_file_structure()
    
    print("🎯 Ready to visualize your SQANTI3 results with UCSC Genome Browser!")

if __name__ == "__main__":
    main()
