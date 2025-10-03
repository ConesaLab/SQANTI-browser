#!/bin/bash
# UCSC Tools Installation Script for SQANTI3 to UCSC Genome Browser Integration
# This script downloads and installs the required UCSC tools

set -e

echo "Installing UCSC tools for SQANTI3 to UCSC Genome Browser integration..."

# Detect OS and architecture
OS=$(uname -s)
ARCH=$(uname -m)

echo "Detected OS: $OS, Architecture: $ARCH"

# Set download URL based on OS
if [[ "$OS" == "Linux" ]]; then
    if [[ "$ARCH" == "x86_64" ]]; then
        BASE_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64"
    elif [[ "$ARCH" == "aarch64" ]]; then
        BASE_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.aarch64"
    else
        echo "Unsupported Linux architecture: $ARCH"
        exit 1
    fi
elif [[ "$OS" == "Darwin" ]]; then
    if [[ "$ARCH" == "x86_64" ]]; then
        BASE_URL="http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64"
    elif [[ "$ARCH" == "arm64" ]]; then
        BASE_URL="http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64"
    else
        echo "Unsupported macOS architecture: $ARCH"
        exit 1
    fi
else
    echo "Unsupported operating system: $OS"
    echo "Please install UCSC tools manually from: http://hgdownload.soe.ucsc.edu/admin/exe/"
    exit 1
fi

# Required tools
TOOLS=("gtfToGenePred" "genePredToBed" "bedToBigBed" "hubCheck", "bigBedInfo")

# Create installation directory
INSTALL_DIR="/usr/local/bin"
if [[ ! -w "$INSTALL_DIR" ]]; then
    echo "Warning: Cannot write to $INSTALL_DIR. Installing to ~/bin instead."
    INSTALL_DIR="$HOME/bin"
    mkdir -p "$INSTALL_DIR"
    echo "Please add $INSTALL_DIR to your PATH by adding this line to your shell profile:"
    echo "export PATH=\"$INSTALL_DIR:\$PATH\""
fi

echo "Installing to: $INSTALL_DIR"

# Download and install each tool
for tool in "${TOOLS[@]}"; do
    echo "Downloading $tool..."
    
    # Download the tool
    if curl -L -o "$tool" "$BASE_URL/$tool"; then
        echo "Downloaded $tool successfully"
        
        # Make executable
        chmod +x "$tool"
        
        # Move to installation directory
        if mv "$tool" "$INSTALL_DIR/"; then
            echo "Installed $tool to $INSTALL_DIR"
        else
            echo "Failed to install $tool to $INSTALL_DIR"
            exit 1
        fi
    else
        echo "Failed to download $tool"
        exit 1
    fi
done

echo ""
echo "Installation completed successfully!"
echo ""

# Verify installation
echo "Verifying installation..."
for tool in "${TOOLS[@]}"; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo "✓ $tool is available"
    else
        echo "✗ $tool is not available"
        echo "  Please ensure $INSTALL_DIR is in your PATH"
    fi
done

echo ""
echo "Next steps:"
echo "1. Install Python dependencies: pip install -r requirements.txt"
echo "2. Run the converter: python sqanti3_to_UCSC.py --help"
echo ""
echo "Hub troubleshooting (optional):"
echo "- Validate a hub locally with UCSC hubCheck:"
echo "  hubCheck -noTracks /path/to/hub.txt"
echo ""
echo "If you installed to ~/bin, remember to add it to your PATH:"
echo "export PATH=\"$INSTALL_DIR:\$PATH\""
