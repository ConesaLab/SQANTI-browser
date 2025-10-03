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
    
    print("SQANTI3 to UCSC Genome Browser Integration - Example Usage")
    print("=" * 60)
    
    # Check if the main script exists
    main_script = "sqanti3_to_UCSC.py"
    if not os.path.exists(main_script):
        print(f"Error: {main_script} not found in current directory")
        print("Please ensure you're running this from the directory containing the integration scripts")
        return False
    
    # Example parameters
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
    if not os.path.exists(example_gtf):
        print(f"Note: Example GTF file not found: {example_gtf}")
        print("This is expected if you haven't downloaded the SQANTI3 example data yet.")
        print()
        print("To download example data:")
        print("  git clone https://github.com/ConesaLab/SQANTI3.git")
        print("  cp -r SQANTI3/example/SQANTI3_QC_output ./example/")
        print()
    
    if not os.path.exists(example_classification):
        print(f"Note: Example classification file not found: {example_classification}")
        print("This is expected if you haven't downloaded the SQANTI3 example data yet.")
        print()
    
    # Show the command that would be run
    print("To run the conversion with your own data, use:")
    print()
    print(f"python {main_script} \\")
    print(f"    --gtf <your_corrected.gtf> \\")
    print(f"    --classification <your_classification.txt> \\")
    print(f"    --output <output_directory> \\")
    print(f"    --genome <genome_assembly>")
    print()
    print("For GitHub integration (recommended for UCSC Genome Browser):")
    print(f"python {main_script} \\")
    print(f"    --gtf <your_corrected.gtf> \\")
    print(f"    --classification <your_classification.txt> \\")
    print(f"    --output <output_directory> \\")
    print(f"    --genome <genome_assembly> \\")
    print(f"    --github-repo <username/repository>")
    print()
    print("Optional enhancements:")
    print("- Enable Trix text index (requires ixIxx in PATH): add --enable-trix")
    print("- Add STAR splice junctions track from SJ.out.tab: --star-sj <path>")
    print("- Use a genome .2bit to compute chrom.sizes: --twobit <genome.2bit>")
    print()
    
    # Show help
    print("For detailed help and all options:")
    print(f"python {main_script} --help")
    print()
    
    # Check if we can run the script
    try:
        result = subprocess.run([sys.executable, main_script, "--help"], 
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
    print("3. Upload the generated hub files to a web server")
    print("4. Add the hub to UCSC Genome Browser")
    print("5. View and filter your transcripts by structural category")
    print()
    print("For detailed instructions, see README.md")
    
    return True

def show_file_structure():
    """Show the expected file structure"""
    print("Expected file structure after conversion:")
    print()
    print("output_directory/")
    print("├── hg38_sqanti3.bb                    # Main bigBed file")
    print("├── hub.txt                             # Hub configuration")
    print("├── genomes.txt                         # Genome mapping")
    print("├── hg38/                               # Genome-specific directory")
    print("│   └── trackDb.txt                    # Track configuration")
    print("└── output_directory_hg38_SQANTI3_Hub.html  # Documentation")
    print()
    print("Note: Hub name includes output directory name for uniqueness!")
    print()
    print("Output directory examples:")
    print("  ./my_project          → my_project_hg38_SQANTI3_Hub")
    print("  ~/Documents/research  → research_hg38_SQANTI3_Hub")
    print("  /home/user/analysis   → analysis_hg38_SQANTI3_Hub")
    print("  ~/my_project          → my_project_hg38_SQANTI3_Hub")
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
