#!/usr/bin/env python3
"""
Comprehensive Test Suite for SQANTI-browser

This script tests all major functionalities of SQANTI-browser to ensure
everything works correctly after code changes.

Tests include:
- Installation checks (UCSC tools, Python deps, file permissions, tool functionality)
- Basic GTF + classification conversion
- Sorting options
- Table generation
- Category tracks
- Custom genome (.2bit) support
- STAR splice junctions
- CAGE peaks
- PolyA peaks
- Reference GTF
- Trix search indexing
- File validation

Usage:
  python test_sqanti_browser.py           # Run full test suite
  python test_sqanti_browser.py --install-only   # Run only installation checks (fast)

Author: Carolina Monzo
Date: January 2026
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


class Colors:
    """ANSI color codes for terminal output"""
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    BOLD = '\033[1m'
    END = '\033[0m'


class SQANTIBrowserTester:
    """Test suite for SQANTI-browser"""
    
    def __init__(self):
        self.project_root = Path(__file__).parent
        self.test_output_dir = self.project_root / "test_output"
        self.passed_tests = 0
        self.failed_tests = 0
        self.test_results = []
        
    def print_header(self, text):
        """Print a formatted header"""
        print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.END}")
        print(f"{Colors.BOLD}{Colors.BLUE}{text:^70}{Colors.END}")
        print(f"{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.END}\n")
        
    def print_test(self, test_name):
        """Print test name"""
        print(f"{Colors.BOLD}Testing: {test_name}{Colors.END}")
        
    def print_pass(self, message):
        """Print success message"""
        print(f"  {Colors.GREEN}✓ PASS{Colors.END}: {message}")
        self.passed_tests += 1
        self.test_results.append({"test": message, "status": "PASS"})
        
    def print_fail(self, message, error=None):
        """Print failure message"""
        print(f"  {Colors.RED}✗ FAIL{Colors.END}: {message}")
        if error:
            print(f"    Error: {error}")
        self.failed_tests += 1
        self.test_results.append({"test": message, "status": "FAIL", "error": str(error) if error else None})
        
    def print_warning(self, message):
        """Print warning message"""
        print(f"  {Colors.YELLOW}⚠ WARNING{Colors.END}: {message}")
        
    def cleanup(self):
        """Clean up test output directory"""
        if self.test_output_dir.exists():
            shutil.rmtree(self.test_output_dir)
        self.test_output_dir.mkdir(exist_ok=True)
        
    def check_dependencies(self):
        """Check that all required tools are available"""
        self.print_header("CHECKING DEPENDENCIES")
        
        tools = [
            'bedToBigBed',
            'gtfToGenePred',
            'genePredToBed',
            'ixIxx',
            'twoBitInfo',
            'hubCheck'
        ]
        
        all_present = True
        for tool in tools:
            result = subprocess.run(['which', tool], capture_output=True, text=True)
            if result.returncode == 0:
                self.print_pass(f"{tool} is installed")
            else:
                self.print_fail(f"{tool} is not installed")
                all_present = False
                
        # Check Python dependencies
        try:
            import pandas
            self.print_pass("pandas is installed")
        except ImportError:
            self.print_fail("pandas is not installed")
            all_present = False
            
        return all_present

    def check_file_permissions(self):
        """Check that we can create and write files"""
        self.print_header("CHECKING FILE PERMISSIONS")
        try:
            temp_dir = tempfile.mkdtemp(prefix="sqanti_test_")
            test_file = os.path.join(temp_dir, "test.txt")
            with open(test_file, 'w') as f:
                f.write("test")
            os.remove(test_file)
            os.rmdir(temp_dir)
            self.print_pass("Can create temporary directories and write files")
            return True
        except Exception as e:
            self.print_fail(f"File permission test failed: {e}")
            return False

    def check_ucsc_tool_functionality(self):
        """Test that UCSC tools respond to help commands"""
        self.print_header("CHECKING UCSC TOOL FUNCTIONALITY")
        tools_to_check = [
            ('gtfToGenePred', ['usage:']),
            ('genePredToBed', ['usage:']),
            ('bedToBigBed', ['usage:']),
            ('ixIxx', ['usage:', 'ixixx']),
        ]
        all_ok = True
        for tool, keywords in tools_to_check:
            try:
                result = subprocess.run(
                    [tool],
                    capture_output=True,
                    text=True,
                    timeout=10,
                )
                combined = (result.stdout or '') + (result.stderr or '').lower()
                if any(kw.lower() in combined for kw in keywords):
                    self.print_pass(f"{tool} responds to help command")
                else:
                    self.print_fail(f"{tool} help command failed")
                    all_ok = False
            except subprocess.TimeoutExpired:
                self.print_fail(f"{tool} timed out")
                all_ok = False
            except Exception as e:
                self.print_fail(f"{tool} failed: {e}")
                all_ok = False
        return all_ok

    def run_install_only_tests(self):
        """Run installation checks only (fast, no conversion). Returns True if all pass."""
        self.print_header("INSTALLATION-ONLY TEST MODE")
        all_ok = True
        if not self.check_dependencies():
            all_ok = False
        if not self.check_file_permissions():
            all_ok = False
        if not self.check_ucsc_tool_functionality():
            all_ok = False
        print("\n" + "=" * 60)
        if all_ok:
            print("All installation checks passed. System is ready.")
            print("Run without --install-only for full conversion tests.")
        else:
            print("Some installation checks failed. Fix the issues above.")
            print("Common solutions:")
            print("  1. Install UCSC tools: bash install_ucsc_tools.sh")
            print("  2. Install Python deps: pip install -r requirements.txt")
        return all_ok
        
    def check_example_files(self):
        """Check that all required example files exist"""
        self.print_header("CHECKING EXAMPLE FILES")
        
        required_files = [
            "example/SQANTI3_QC_output/example_corrected.gtf",
            "example/SQANTI3_QC_output/example_classification.txt",
            "example/SQANTI3_QC_output/exampleSJ.out.tab",
            "example/SQANTI3_QC_output/chr22.human.refTSS_v3.1.hg38.bed",
            "example/SQANTI3_QC_output/polyApeaks.atlas.GRCh38.bed",
            "example/SQANTI3_QC_output/gencode.v38.basic_chr22.gtf",
            "example/SQANTI3_QC_output/chrom.sizes",
            "example/SQANTI3_QC_custom_genome/sirv_classification.txt",
            "example/SQANTI3_QC_custom_genome/sirv_corrected.gtf",
            "example/SQANTI3_QC_custom_genome/SIRVS.2bit",
        ]
        
        all_present = True
        for file_path in required_files:
            full_path = self.project_root / file_path
            if full_path.exists():
                self.print_pass(f"{file_path} exists")
            else:
                self.print_fail(f"{file_path} not found")
                all_present = False
                
        return all_present
        
    def run_sqanti_browser(self, test_name, output_subdir, args):
        """Run sqanti3_to_UCSC.py with given arguments"""
        output_dir = self.test_output_dir / output_subdir
        
        cmd = [
            sys.executable,
            str(self.project_root / "sqanti3_to_UCSC.py"),
            "--output", str(output_dir)
        ] + args
        
        print(f"  Running: {' '.join([str(c) for c in cmd])}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode == 0:
                self.print_pass(f"{test_name} completed successfully")
                return True, output_dir
            else:
                self.print_fail(f"{test_name} failed", result.stderr)
                return False, output_dir
                
        except subprocess.TimeoutExpired:
            self.print_fail(f"{test_name} timed out after 5 minutes")
            return False, output_dir
        except Exception as e:
            self.print_fail(f"{test_name} raised exception", str(e))
            return False, output_dir
            
    def validate_hub_structure(self, output_dir, genome, required_files):
        """Validate that hub has correct structure"""
        
        # Check root files
        root_files = ["hub.txt", "genomes.txt", "README.md"]
        for file in root_files:
            path = output_dir / file
            if path.exists():
                self.print_pass(f"Found {file}")
            else:
                self.print_fail(f"Missing {file}")
                
        # Check genome directory
        genome_dir = output_dir / genome
        if not genome_dir.exists():
            self.print_fail(f"Genome directory {genome} not found")
            return False
            
        self.print_pass(f"Found genome directory {genome}/")
        
        # Check required files in genome directory
        for file in required_files:
            path = genome_dir / file
            if path.exists() and path.stat().st_size > 0:
                self.print_pass(f"Found {genome}/{file} ({path.stat().st_size} bytes)")
            else:
                self.print_fail(f"Missing or empty {genome}/{file}")
                
        return True
        
    def validate_hub_with_hubcheck(self, output_dir):
        """Validate hub using UCSC hubCheck"""
        hub_txt = output_dir / "hub.txt"
        
        if not hub_txt.exists():
            self.print_fail("hub.txt not found, cannot run hubCheck")
            return False
            
        try:
            result = subprocess.run(
                ['hubCheck', str(hub_txt)],
                capture_output=True,
                text=True,
                timeout=60
            )
            
            if result.returncode == 0:
                self.print_pass("hubCheck validation passed")
                return True
            else:
                self.print_fail("hubCheck validation failed", result.stderr)
                return False
                
        except subprocess.TimeoutExpired:
            self.print_fail("hubCheck timed out")
            return False
        except Exception as e:
            self.print_fail("hubCheck raised exception", str(e))
            return False
            
    def test_basic_conversion(self):
        """Test 1: Basic GTF + classification conversion"""
        self.print_header("TEST 1: Basic Conversion")
        self.print_test("Basic GTF + classification conversion")
        
        success, output_dir = self.run_sqanti_browser(
            "Basic conversion",
            "test1_basic",
            [
                "--gtf", "example/SQANTI3_QC_output/example_corrected.gtf",
                "--classification", "example/SQANTI3_QC_output/example_classification.txt",
                "--genome", "hg38",
                "--chrom-sizes", "example/SQANTI3_QC_output/chrom.sizes"
            ]
        )
        
        if success:
            required_files = [
                "hg38_sqanti3.bb",
                "hg38_sqanti3_track.html",
                "trackDb.txt",
                "groups.txt",
                "trix.ix",
                "trix.ixx"
            ]
            self.validate_hub_structure(output_dir, "hg38", required_files)
            
        return success
        
    def test_with_category_tracks(self):
        """Test 2: Conversion with category-specific tracks"""
        self.print_header("TEST 2: Category-Specific Tracks")
        self.print_test("Conversion with per-category tracks")
        
        success, output_dir = self.run_sqanti_browser(
            "Category tracks",
            "test2_categories",
            [
                "--gtf", "example/SQANTI3_QC_output/example_corrected.gtf",
                "--classification", "example/SQANTI3_QC_output/example_classification.txt",
                "--genome", "hg38",
                "--chrom-sizes", "example/SQANTI3_QC_output/chrom.sizes"
            ]
        )
        
        if success:
            genome_dir = output_dir / "hg38"
            
            # Check for category-specific bigBed files
            categories = [
                "full-splice_match",
                "incomplete-splice_match",
                "novel_in_catalog",
                "novel_not_in_catalog"
            ]
            
            for category in categories:
                bb_file = genome_dir / f"hg38_sqanti3_{category}.bb"
                html_file = genome_dir / f"hg38_sqanti3_{category}.html"
                
                if bb_file.exists():
                    self.print_pass(f"Found category bigBed: {category}")
                else:
                    self.print_fail(f"Missing category bigBed: {category}")
                    
                if html_file.exists():
                    self.print_pass(f"Found category HTML: {category}")
                else:
                    self.print_fail(f"Missing category HTML: {category}")
                    
        return success
        
    def test_with_tables(self):
        """Test 3: Conversion with HTML tables"""
        self.print_header("TEST 3: HTML Tables")
        self.print_test("Conversion with interactive HTML tables")
        
        success, output_dir = self.run_sqanti_browser(
            "HTML tables",
            "test3_tables",
            [
                "--gtf", "example/SQANTI3_QC_output/example_corrected.gtf",
                "--classification", "example/SQANTI3_QC_output/example_classification.txt",
                "--genome", "hg38",
                "--chrom-sizes", "example/SQANTI3_QC_output/chrom.sizes",
                "--tables"
            ]
        )
        
        if success:
            tables_dir = output_dir / "table_reports"
            
            if tables_dir.exists():
                self.print_pass("table_reports directory created")
                
                # Check for at least some HTML tables
                html_files = list(tables_dir.glob("*.html"))
                if html_files:
                    self.print_pass(f"Found {len(html_files)} HTML table files")
                else:
                    self.print_fail("No HTML table files found")
            else:
                self.print_fail("table_reports directory not created")
                
        return success
        
    def test_sorting_options(self):
        """Test 4: Different sorting options"""
        self.print_header("TEST 4: Sorting Options")
        
        sort_options = ["none", "length", "iso_exp"]
        all_passed = True
        
        for sort_by in sort_options:
            self.print_test(f"Sorting by: {sort_by}")
            
            success, output_dir = self.run_sqanti_browser(
                f"Sort by {sort_by}",
                f"test4_sort_{sort_by}",
                [
                    "--gtf", "example/SQANTI3_QC_output/example_corrected.gtf",
                    "--classification", "example/SQANTI3_QC_output/example_classification.txt",
                    "--genome", "hg38",
                    "--chrom-sizes", "example/SQANTI3_QC_output/chrom.sizes",
                    "--sort-by", sort_by,
                    "--no-category-tracks"  # Speed up test
                ]
            )
            
            if not success:
                all_passed = False
                
        return all_passed
        
    def test_custom_genome_twobit(self):
        """Test 5: Custom genome with .2bit file"""
        self.print_header("TEST 5: Custom Genome (.2bit)")
        self.print_test("Conversion with .2bit genome file")
        
        success, output_dir = self.run_sqanti_browser(
            "2bit genome",
            "test5_twobit",
            [
                "--gtf", "example/SQANTI3_QC_custom_genome/sirv_corrected.gtf",
                "--classification", "example/SQANTI3_QC_custom_genome/sirv_classification.txt",
                "--genome", "SIRVS",
                "--twobit", "example/SQANTI3_QC_custom_genome/SIRVS.2bit",
                "--no-category-tracks"  # Speed up test
            ]
        )
        
        if success:
            # Check that chrom.sizes was auto-generated
            genome_dir = output_dir / "SIRVS"
            
            # Check for 2bit file in genome directory
            twobit_file = genome_dir / "SIRVS.2bit"
            if twobit_file.exists():
                self.print_pass("2bit file copied to genome directory")
            else:
                self.print_fail("2bit file not found in genome directory")
                
            # Check genomes.txt has twoBitPath
            genomes_txt = output_dir / "genomes.txt"
            if genomes_txt.exists():
                with open(genomes_txt) as f:
                    content = f.read()
                    if "twoBitPath" in content:
                        self.print_pass("genomes.txt contains twoBitPath")
                    else:
                        self.print_fail("genomes.txt missing twoBitPath")
                        
                    if "htmlPath" in content:
                        self.print_pass("genomes.txt contains htmlPath")
                    else:
                        self.print_fail("genomes.txt missing htmlPath")
                        
            # Check for genome HTML description
            genome_html = genome_dir / "SIRVS_genome.html"
            if genome_html.exists():
                self.print_pass("Genome description HTML created")
            else:
                self.print_fail("Genome description HTML not found")
                
        return success
        
    def test_validation_tracks(self):
        """Test 6: All validation tracks (STAR-SJ, CAGE, PolyA, RefGTF)"""
        self.print_header("TEST 6: Validation Tracks")
        self.print_test("Conversion with all validation tracks")
        
        success, output_dir = self.run_sqanti_browser(
            "Validation tracks",
            "test6_validation",
            [
                "--gtf", "example/SQANTI3_QC_output/example_corrected.gtf",
                "--classification", "example/SQANTI3_QC_output/example_classification.txt",
                "--genome", "hg38",
                "--chrom-sizes", "example/SQANTI3_QC_output/chrom.sizes",
                "--star-sj", "example/SQANTI3_QC_output/exampleSJ.out.tab",
                "--CAGE-peak", "example/SQANTI3_QC_output/chr22.human.refTSS_v3.1.hg38.bed",
                "--polyA-peak", "example/SQANTI3_QC_output/polyApeaks.atlas.GRCh38.bed",
                "--refGTF", "example/SQANTI3_QC_output/gencode.v38.basic_chr22.gtf",
                "--no-category-tracks"  # Speed up test
            ]
        )
        
        if success:
            genome_dir = output_dir / "hg38"
            
            validation_files = {
                "hg38_star_sj.bb": "STAR splice junctions",
                "hg38_cage_peaks.bb": "CAGE peaks",
                "hg38_polya_peaks.bb": "PolyA peaks",
                "hg38_reference.bb": "Reference GTF"
            }
            
            for file, description in validation_files.items():
                path = genome_dir / file
                if path.exists() and path.stat().st_size > 0:
                    self.print_pass(f"Found {description}: {file}")
                else:
                    self.print_fail(f"Missing {description}: {file}")
                    
            # Check groups.txt has validation group
            groups_txt = genome_dir / "groups.txt"
            if groups_txt.exists():
                with open(groups_txt) as f:
                    content = f.read()
                    if "validation" in content.lower():
                        self.print_pass("Validation group defined in groups.txt")
                    else:
                        self.print_fail("Validation group not found in groups.txt")
                        
        return success
        
    def test_no_category_tracks(self):
        """Test 7: Main track only (no category tracks)"""
        self.print_header("TEST 7: Main Track Only")
        self.print_test("Conversion without category-specific tracks")
        
        success, output_dir = self.run_sqanti_browser(
            "No category tracks",
            "test7_no_categories",
            [
                "--gtf", "example/SQANTI3_QC_output/example_corrected.gtf",
                "--classification", "example/SQANTI3_QC_output/example_classification.txt",
                "--genome", "hg38",
                "--chrom-sizes", "example/SQANTI3_QC_output/chrom.sizes",
                "--no-category-tracks"
            ]
        )
        
        if success:
            genome_dir = output_dir / "hg38"
            
            # Should have main track
            main_bb = genome_dir / "hg38_sqanti3.bb"
            if main_bb.exists():
                self.print_pass("Main bigBed file created")
            else:
                self.print_fail("Main bigBed file not found")
                
            # Should NOT have category tracks
            category_bbs = list(genome_dir.glob("hg38_sqanti3_*.bb"))
            # Filter out validation tracks
            category_bbs = [f for f in category_bbs if not any(x in f.name for x in ['star', 'cage', 'polya', 'reference'])]
            
            if not category_bbs:
                self.print_pass("No category-specific tracks created (as expected)")
            else:
                self.print_fail(f"Unexpected category tracks found: {[f.name for f in category_bbs]}")
                
        return success
        
    def test_dry_run(self):
        """Test 8: Dry run mode"""
        self.print_header("TEST 8: Dry Run Mode")
        self.print_test("Dry run (no bigBed conversion)")
        
        success, output_dir = self.run_sqanti_browser(
            "Dry run",
            "test8_dryrun",
            [
                "--gtf", "example/SQANTI3_QC_output/example_corrected.gtf",
                "--classification", "example/SQANTI3_QC_output/example_classification.txt",
                "--genome", "hg38",
                "--chrom-sizes", "example/SQANTI3_QC_output/chrom.sizes",
                "--dry-run",
                "--no-category-tracks"
            ]
        )
        
        if success:
            # In dry run, should NOT have bigBed files or hub files
            genome_dir = output_dir / "hg38"
            hub_txt = output_dir / "hub.txt"
            
            if not hub_txt.exists():
                self.print_pass("hub.txt not created in dry run (as expected)")
            else:
                self.print_warning("hub.txt was created in dry run mode")
                
        return success
        
    def test_hubcheck_validation(self):
        """Test 9: Validate hub with hubCheck"""
        self.print_header("TEST 9: Hub Validation (hubCheck)")
        self.print_test("Running hubCheck on generated hub")
        
        # Use a previously generated hub (test 1)
        test1_output = self.test_output_dir / "test1_basic"
        
        if not test1_output.exists():
            self.print_warning("Test 1 output not found, running basic conversion first")
            self.test_basic_conversion()
            
        if test1_output.exists():
            return self.validate_hub_with_hubcheck(test1_output)
        else:
            self.print_fail("Could not validate hub - no output available")
            return False
            
    def run_all_tests(self):
        """Run all tests"""
        start_time = time.time()
        
        print(f"\n{Colors.BOLD}{Colors.BLUE}")
        print("╔════════════════════════════════════════════════════════════════════╗")
        print("║                                                                    ║")
        print("║              SQANTI-browser Comprehensive Test Suite              ║")
        print("║                                                                    ║")
        print("╚════════════════════════════════════════════════════════════════════╝")
        print(Colors.END)
        
        # Pre-flight checks
        if not self.check_dependencies():
            print(f"\n{Colors.RED}{Colors.BOLD}ERROR: Missing required dependencies{Colors.END}")
            print("Please run: bash install_ucsc_tools.sh")
            print("And: pip install -r requirements.txt")
            return False
            
        if not self.check_example_files():
            print(f"\n{Colors.RED}{Colors.BOLD}ERROR: Missing required example files{Colors.END}")
            print("Please ensure all example files are present")
            return False
            
        # Clean up old test outputs
        self.print_header("PREPARING TEST ENVIRONMENT")
        print("Cleaning up old test outputs...")
        self.cleanup()
        self.print_pass("Test environment ready")
        
        # Run all tests
        tests = [
            ("Basic Conversion", self.test_basic_conversion),
            ("Category Tracks", self.test_with_category_tracks),
            ("HTML Tables", self.test_with_tables),
            ("Sorting Options", self.test_sorting_options),
            ("Custom Genome (.2bit)", self.test_custom_genome_twobit),
            ("Validation Tracks", self.test_validation_tracks),
            ("Main Track Only", self.test_no_category_tracks),
            ("Dry Run Mode", self.test_dry_run),
            ("Hub Validation", self.test_hubcheck_validation),
        ]
        
        for test_name, test_func in tests:
            try:
                test_func()
            except Exception as e:
                self.print_fail(f"{test_name} raised unexpected exception", str(e))
                
        # Print summary
        elapsed_time = time.time() - start_time
        self.print_summary(elapsed_time)
        
        return self.failed_tests == 0
        
    def print_summary(self, elapsed_time):
        """Print test summary"""
        total_tests = self.passed_tests + self.failed_tests
        pass_rate = (self.passed_tests / total_tests * 100) if total_tests > 0 else 0
        
        print(f"\n{Colors.BOLD}{Colors.BLUE}")
        print("╔════════════════════════════════════════════════════════════════════╗")
        print("║                          TEST SUMMARY                              ║")
        print("╚════════════════════════════════════════════════════════════════════╝")
        print(Colors.END)
        
        print(f"\n{Colors.BOLD}Results:{Colors.END}")
        print(f"  Total Tests: {total_tests}")
        print(f"  {Colors.GREEN}Passed: {self.passed_tests}{Colors.END}")
        print(f"  {Colors.RED}Failed: {self.failed_tests}{Colors.END}")
        print(f"  Pass Rate: {pass_rate:.1f}%")
        print(f"  Time Elapsed: {elapsed_time:.1f} seconds")
        
        if self.failed_tests == 0:
            print(f"\n{Colors.GREEN}{Colors.BOLD}✓ ALL TESTS PASSED!{Colors.END}")
            print(f"{Colors.GREEN}SQANTI-browser is working correctly.{Colors.END}\n")
        else:
            print(f"\n{Colors.RED}{Colors.BOLD}✗ SOME TESTS FAILED{Colors.END}")
            print(f"{Colors.RED}Please review the errors above.{Colors.END}\n")
            
        # Save results to JSON
        results_file = self.test_output_dir / "test_results.json"
        with open(results_file, 'w') as f:
            json.dump({
                "total_tests": total_tests,
                "passed": self.passed_tests,
                "failed": self.failed_tests,
                "pass_rate": pass_rate,
                "elapsed_time": elapsed_time,
                "results": self.test_results
            }, f, indent=2)
            
        print(f"Detailed results saved to: {results_file}")
        print(f"Test outputs saved to: {self.test_output_dir}/\n")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="SQANTI-browser test suite. Run full tests or installation-only checks."
    )
    parser.add_argument(
        "--install-only",
        action="store_true",
        help="Run only installation checks (UCSC tools, Python deps, file permissions). Fast.",
    )
    args = parser.parse_args()

    tester = SQANTIBrowserTester()
    if args.install_only:
        success = tester.run_install_only_tests()
    else:
        success = tester.run_all_tests()

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
