# example_output_hg38_SQANTI3_Hub

This hub displays SQANTI3 transcriptome analysis results for the hg38 genome assembly.

## 🚀 Usage Instructions

1. Upload all files in this directory to a web-accessible location (e.g., GitHub).
2. In the UCSC Genome Browser, go to **My Data → Track Hubs**.
3. Enter the URL to your `hub.txt` file.
4. Select the appropriate genome assembly (hg38).
5. The SQANTI3 tracks will appear in your track list.

## 🔍 Advanced Filtering

This hub uses the bigBed 12+44 format with native UCSC filters. Right-click on the track and select "Configure" or "Filter" to access controls.

## 🔎 Trix Search

Use the search box to find isoforms by attribute. Search terms use underscores (e.g., `structural_category_full_splice_match`).

## 🎨 Color Legend

### Standard Colors

- **Full-splice Match (FSM):** #6BAED6 (RGB 107, 174, 214)
- **Incomplete-splice Match (ISM):** #FC8D59 (RGB 252, 141, 89)
- **Novel In Catalog (NIC):** #78C679 (RGB 120, 198, 121)
- **Novel Not In Catalog (NNC):** #D62F4B (RGB 214, 47, 75)
- **Genic:** #969696 (RGB 150, 150, 150)
- **Antisense:** #66C2A4 (RGB 102, 194, 164)
- **Fusion:** #DAA520 (RGB 218, 165, 32)
- **Intergenic:** #E9967A (RGB 233, 150, 122)
- **Genic Intron:** #41B6C4 (RGB 65, 182, 196)

### Highlight Colors (top FL isoform per group)

- **Full-splice Match (FSM):** #456F89 (RGB 69, 111, 137)
- **Incomplete-splice Match (ISM):** #CA7147 (RGB 202, 113, 71)
- **Novel In Catalog (NIC):** #4D7E4E (RGB 77, 126, 78)
- **Novel Not In Catalog (NNC):** #98071F (RGB 152, 7, 31)
- **Genic:** #606060 (RGB 96, 96, 96)
- **Antisense:** #427C69 (RGB 66, 124, 105)
- **Fusion:** #8B6A15 (RGB 139, 106, 21)
- **Intergenic:** #95604E (RGB 149, 96, 78)
- **Genic Intron:** #2A757E (RGB 42, 117, 126)
