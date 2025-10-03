#!/usr/bin/env python3
"""
Test script to validate SQANTI3 to UCSC Genome Browser integration setup, it checks that all required tools and dependencies are properly installed.

Author: Carolina Monzo
"""

import subprocess
import sys
import shutil
import importlib

def test_ucsc_tools():
    """Test if UCSC tools are available"""
    print("Testing UCSC tools...")
    
    required_tools = ['gtfToGenePred', 'genePredToBed', 'bedToBigBed']
    missing_tools = []
    
    for tool in required_tools:
        if shutil.which(tool) is None:
            missing_tools.append(tool)
            print(f"  ✗ {tool} - NOT FOUND")
        else:
            print(f"  ✓ {tool} - Found at {shutil.which(tool)}")
    
    if missing_tools:
        print(f"\nMissing tools: {', '.join(missing_tools)}")
        print("Please run: bash install_ucsc_tools.sh")
        return False
    
    print("  All UCSC tools are available!")
    return True

def test_python_dependencies():
    """Test if Python dependencies are available"""
    print("\nTesting Python dependencies...")
    
    required_modules = ['pandas', 'pathlib']
    missing_modules = []
    
    for module in required_modules:
        try:
            importlib.import_module(module)
            print(f"  ✓ {module} - Available")
        except ImportError:
            missing_modules.append(module)
            print(f"  ✗ {module} - NOT FOUND")
    
    if missing_modules:
        print(f"\nMissing Python modules: {', '.join(missing_modules)}")
        print("Please run: pip install -r requirements.txt")
        return False
    
    print("  All Python dependencies are available!")
    return True

def test_file_permissions():
    """Test if we can create and write files"""
    print("\nTesting file permissions...")
    
    try:
        import tempfile
        import os
        
        # Test temp directory creation
        temp_dir = tempfile.mkdtemp(prefix="test_")
        print(f"  ✓ Can create temporary directories: {temp_dir}")
        
        # Test file writing
        test_file = os.path.join(temp_dir, "test.txt")
        with open(test_file, 'w') as f:
            f.write("test")
        print("  ✓ Can write files")
        
        # Cleanup
        os.remove(test_file)
        os.rmdir(temp_dir)
        print("  ✓ Can delete files and directories")
        
        return True
        
    except Exception as e:
        print(f"  ✗ File permission test failed: {e}")
        return False

def test_ucsc_tool_functionality():
    """Test if UCSC tools can run basic commands"""
    print("\nTesting UCSC tool functionality...")
    
    try:
        # Test gtfToGenePred help
        result = subprocess.run(['gtfToGenePred'], 
                              capture_output=True, text=True, timeout=10)
        if "usage:" in result.stdout or "usage:" in result.stderr:
            print("  ✓ gtfToGenePred responds to help command")
        else:
            print("  ✗ gtfToGenePred help command failed")
            return False
        
        # Test genePredToBed help
        result = subprocess.run(['genePredToBed'], 
                              capture_output=True, text=True, timeout=10)
        if "usage:" in result.stdout or "usage:" in result.stderr:
            print("  ✓ genePredToBed responds to help command")
        else:
            print("  ✗ genePredToBed help command failed")
            return False
        
        # Test bedToBigBed help
        result = subprocess.run(['bedToBigBed'], 
                              capture_output=True, text=True, timeout=10)
        if "usage:" in result.stdout or "usage:" in result.stderr:
            print("  ✓ bedToBigBed responds to help command")
        else:
            print("  ✗ bedToBigBed help command failed")
            return False
        
        return True
        
    except subprocess.TimeoutExpired:
        print("  ✗ UCSC tool help commands timed out")
        return False
    except Exception as e:
        print(f"  ✗ UCSC tool functionality test failed: {e}")
        return False

def main():
    """Run all tests"""
    print("SQANTI3 to UCSC Genome Browser Integration - Installation Test")
    print("=" * 60)
    
    all_tests_passed = True
    
    # Run tests
    if not test_ucsc_tools():
        all_tests_passed = False
    
    if not test_python_dependencies():
        all_tests_passed = False
    
    if not test_file_permissions():
        all_tests_passed = False
    
    if not test_ucsc_tool_functionality():
        all_tests_passed = False
    
    # Summary
    print("\n" + "=" * 60)
    if all_tests_passed:
        print("🎉 All tests passed! Your system is ready for SQANTI3 to UCSC integration.")
        print("\nYou can now run:")
        print("  python sqanti3_to_UCSC.py --help")
    else:
        print("❌ Some tests failed. Please fix the issues above before proceeding.")
        print("\nCommon solutions:")
        print("  1. Install UCSC tools: bash install_ucsc_tools.sh")
        print("  2. Install Python dependencies: pip install -r requirements.txt")
        print("  3. Check file permissions and PATH settings")
    
    print("\nFor detailed help, see README.md")
    
    return 0 if all_tests_passed else 1

if __name__ == "__main__":
    sys.exit(main())
