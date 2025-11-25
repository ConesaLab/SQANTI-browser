# example_output_hg38_SQANTI3_Hub

This hub displays SQANTI3 transcriptome analysis results for the hg38 genome assembly.

## 🔍 Advanced Filtering

**This hub uses the bigBed 12+44 format with native UCSC filters.**

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
4. Select the appropriate genome assembly (hg38).
5. The SQANTI3 tracks will appear in your track list.

## 📧 Contact Information

**Note:** The `hub.txt` file includes a placeholder or GitHub-specific email address. If you need to use a different email:
- Edit the `hub.txt` file.
- Change the `email` line to your preferred email address.
- Re-upload the updated file.
