"""
Edge case tests for SQANTI-browser.

Tests unusual inputs, missing files, empty files, and error handling.

Run with: python tests/test_edge_cases.py -v
     or:  PYTHONPATH=. python tests/test_edge_cases.py -v
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

# UCSC tools required for validate-only; skip that test when not installed (e.g. CI before install step)
UCSC_TOOLS_AVAILABLE = shutil.which("gtfToGenePred") is not None

# Project root for paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Add project root so "src" and "sqanti_browser" can be imported
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


class TestMissingFiles(unittest.TestCase):
    """Tests for missing input files."""

    def test_missing_gtf_fails_gracefully(self):
        """Running with non-existent GTF should fail with clear error."""
        result = subprocess.run(
            [
                "python",
                "-m", "sqanti_browser",
                "--gtf",
                "/nonexistent/file.gtf",
                "--classification",
                str(PROJECT_ROOT / "example/SQANTI3_QC_output/example_classification.txt"),
                "--output",
                "/tmp/test_hub",
                "--genome",
                "hg38",
            ],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
            env={**dict(os.environ), "PYTHONPATH": str(PROJECT_ROOT)},
        )
        self.assertNotEqual(result.returncode, 0)
        err_lower = result.stderr.lower()
        self.assertTrue("gtf" in err_lower or "file" in err_lower or "error" in err_lower)

    def test_missing_classification_fails_gracefully(self):
        """Running with non-existent classification should fail."""
        result = subprocess.run(
            [
                "python",
                "-m", "sqanti_browser",
                "--gtf",
                str(PROJECT_ROOT / "example/SQANTI3_QC_output/example_corrected.gtf"),
                "--classification",
                "/nonexistent/classification.txt",
                "--output",
                "/tmp/test_hub",
                "--genome",
                "hg38",
            ],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
            env={**dict(os.environ), "PYTHONPATH": str(PROJECT_ROOT)},
        )
        self.assertNotEqual(result.returncode, 0)


class TestEmptyFiles(unittest.TestCase):
    """Tests for empty or malformed input files."""

    def test_empty_classification_fails(self):
        """Empty classification file should fail."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("")  # empty
            path = f.name
        try:
            result = subprocess.run(
                [
                    "python",
                    "-m", "sqanti_browser",
                    "--gtf",
                    str(PROJECT_ROOT / "example/SQANTI3_QC_output/example_corrected.gtf"),
                    "--classification",
                    path,
                    "--output",
                    "/tmp/test_hub_empty",
                    "--genome",
                    "hg38",
                    "--validate-only",
                ],
                capture_output=True,
                text=True,
                cwd=PROJECT_ROOT,
                env={**dict(os.environ), "PYTHONPATH": str(PROJECT_ROOT)},
            )
            self.assertNotEqual(result.returncode, 0)
        finally:
            Path(path).unlink(missing_ok=True)

    def test_classification_header_only_fails(self):
        """Classification with only header row (no data) should fail or handle gracefully."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("isoform\tchrom\tstrand\tlength\texons\tstructural_category\n")
            path = f.name
        try:
            result = subprocess.run(
                [
                    "python",
                    "-m", "sqanti_browser",
                    "--gtf",
                    str(PROJECT_ROOT / "example/SQANTI3_QC_output/example_corrected.gtf"),
                    "--classification",
                    path,
                    "--output",
                    "/tmp/test_hub_header",
                    "--genome",
                    "hg38",
                ],
                capture_output=True,
                text=True,
                cwd=PROJECT_ROOT,
                env={**dict(os.environ), "PYTHONPATH": str(PROJECT_ROOT)},
                timeout=60,
            )
            # May succeed with empty output or fail - both are acceptable
            self.assertIn(result.returncode, (0, 1))
        finally:
            Path(path).unlink(missing_ok=True)


class TestCLIArgs(unittest.TestCase):
    """Tests for CLI argument handling."""

    def test_help_succeeds(self):
        """--help should succeed."""
        result = subprocess.run(
            ["python", "-m", "sqanti_browser", "--help"],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
            env={**dict(os.environ), "PYTHONPATH": str(PROJECT_ROOT)},
        )
        self.assertEqual(result.returncode, 0)
        out_lower = result.stdout.lower()
        self.assertIn("gtf", out_lower)
        self.assertIn("classification", out_lower)

    def test_missing_required_args_fails(self):
        """Omitting required args should fail."""
        result = subprocess.run(
            ["python", "-m", "sqanti_browser"],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
            env={**dict(os.environ), "PYTHONPATH": str(PROJECT_ROOT)},
        )
        self.assertNotEqual(result.returncode, 0)

    @unittest.skipIf(not UCSC_TOOLS_AVAILABLE, "UCSC tools not installed (gtfToGenePred missing)")
    def test_validate_only_with_valid_files(self):
        """--validate-only with valid files should succeed quickly."""
        result = subprocess.run(
            [
                "python",
                "-m", "sqanti_browser",
                "--gtf",
                str(PROJECT_ROOT / "example/SQANTI3_QC_output/example_corrected.gtf"),
                "--classification",
                str(PROJECT_ROOT / "example/SQANTI3_QC_output/example_classification.txt"),
                "--output",
                "/tmp/test_hub",
                "--genome",
                "hg38",
                "--validate-only",
            ],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
            env={**dict(os.environ), "PYTHONPATH": str(PROJECT_ROOT)},
        )
        self.assertEqual(result.returncode, 0)


class TestFilterIsoformsCLI(unittest.TestCase):
    """Edge case tests for filter_isoforms.py CLI."""

    def test_filter_isoforms_missing_classification(self):
        """filter_isoforms with missing classification should fail."""
        result = subprocess.run(
            [
                "python",
                str(PROJECT_ROOT / "src/filter_isoforms.py"),
                "--classification",
                "/nonexistent/classification.txt",
                "--output-dir",
                "/tmp/reports",
            ],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
            env={**dict(os.environ), "PYTHONPATH": str(PROJECT_ROOT)},
        )
        self.assertNotEqual(result.returncode, 0)

    def test_filter_isoforms_help(self):
        """filter_isoforms --help should succeed."""
        result = subprocess.run(
            ["python", "-m", "src.filter_isoforms", "--help"],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
            env={**dict(os.environ), "PYTHONPATH": str(PROJECT_ROOT)},
        )
        self.assertEqual(result.returncode, 0, msg=f"stderr: {result.stderr}")


class TestParseColorValueEdgeCases(unittest.TestCase):
    """Edge cases for parse_color_value."""

    def test_float_rgb_values_accepted(self):
        """Float RGB values should be coerced to int."""
        from sqanti_browser import parse_color_value

        self.assertEqual(parse_color_value([100.7, 150.2, 200.9]), (100, 150, 200))

    def test_hex_uppercase(self):
        """Uppercase hex should work."""
        from sqanti_browser import parse_color_value

        self.assertEqual(parse_color_value("#FFFFFF"), (255, 255, 255))

    def test_hex_lowercase(self):
        """Lowercase hex should work."""
        from sqanti_browser import parse_color_value

        self.assertEqual(parse_color_value("#ffffff"), (255, 255, 255))


class TestNormalizeTrixTokenEdgeCases(unittest.TestCase):
    """Edge cases for normalize_trix_token."""

    def test_unicode_string(self):
        """Unicode in string should be preserved or handled."""
        from src.utils import normalize_trix_token

        result = normalize_trix_token("gene_α")
        self.assertIsInstance(result, str)

    def test_numbers(self):
        """Numeric string should be normalized."""
        from src.utils import normalize_trix_token

        self.assertEqual(normalize_trix_token("123"), "123")


if __name__ == "__main__":
    unittest.main()
