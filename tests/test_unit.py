"""
Unit tests for SQANTI-browser modules.

Run with: python tests/test_unit.py -v
     or:  PYTHONPATH=. python tests/test_unit.py -v
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

# Add project root so "src" and "sqanti_browser" can be imported
_project_root = Path(__file__).resolve().parent.parent
if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))

import pandas as pd

# ---------------------------------------------------------------------------
# Utils
# ---------------------------------------------------------------------------


class TestNormalizeTrixToken(unittest.TestCase):
    """Tests for src.utils.normalize_trix_token."""

    def test_none_returns_empty(self):
        from src.utils import normalize_trix_token

        self.assertEqual(normalize_trix_token(None), "")

    def test_empty_string_returns_empty(self):
        from src.utils import normalize_trix_token

        self.assertEqual(normalize_trix_token(""), "")
        self.assertEqual(normalize_trix_token("   "), "")

    def test_plus_to_plus(self):
        from src.utils import normalize_trix_token

        self.assertEqual(normalize_trix_token("+"), "plus")

    def test_minus_to_minus(self):
        from src.utils import normalize_trix_token

        self.assertEqual(normalize_trix_token("-"), "minus")

    def test_hyphens_to_underscores(self):
        from src.utils import normalize_trix_token

        self.assertEqual(normalize_trix_token("full-splice_match"), "full_splice_match")

    def test_spaces_to_underscores(self):
        from src.utils import normalize_trix_token

        self.assertEqual(normalize_trix_token("full splice match"), "full_splice_match")

    def test_mixed_special_chars(self):
        from src.utils import normalize_trix_token

        self.assertEqual(normalize_trix_token("structural_category: full-splice_match"), "structural_category_full_splice_match")

    def test_already_normalized(self):
        from src.utils import normalize_trix_token

        self.assertEqual(normalize_trix_token("structural_category_full_splice_match"), "structural_category_full_splice_match")


class TestDarkenColor(unittest.TestCase):
    """Tests for src.utils.darken_color."""

    def test_basic_darken(self):
        from src.utils import darken_color

        self.assertEqual(darken_color((100, 100, 100), factor=0.5), (50, 50, 50))

    def test_default_factor(self):
        from src.utils import darken_color

        r, g, b = darken_color((100, 100, 100))
        self.assertEqual((r, g, b), (64, 64, 64))

    def test_zero_factor(self):
        from src.utils import darken_color

        self.assertEqual(darken_color((255, 255, 255), factor=0), (0, 0, 0))


# ---------------------------------------------------------------------------
# Parse color value (sqanti_browser)
# ---------------------------------------------------------------------------


class TestParseColorValue(unittest.TestCase):
    """Tests for parse_color_value in sqanti_browser."""

    def test_rgb_list(self):
        from sqanti_browser import parse_color_value

        self.assertEqual(parse_color_value([100, 150, 200]), (100, 150, 200))

    def test_rgb_tuple(self):
        from sqanti_browser import parse_color_value

        self.assertEqual(parse_color_value((100, 150, 200)), (100, 150, 200))

    def test_hex_with_hash(self):
        from sqanti_browser import parse_color_value

        self.assertEqual(parse_color_value("#64FF00"), (100, 255, 0))

    def test_hex_without_hash(self):
        from sqanti_browser import parse_color_value

        self.assertEqual(parse_color_value("0000FF"), (0, 0, 255))

    def test_rgb_out_of_range_raises(self):
        from sqanti_browser import parse_color_value

        with self.assertRaises(ValueError):
            parse_color_value([300, 0, 0])

    def test_rgb_wrong_length_raises(self):
        from sqanti_browser import parse_color_value

        with self.assertRaises(ValueError):
            parse_color_value([100, 200])

    def test_hex_wrong_length_raises(self):
        from sqanti_browser import parse_color_value

        with self.assertRaises(ValueError):
            parse_color_value("#FFF")

    def test_invalid_type_raises(self):
        from sqanti_browser import parse_color_value

        with self.assertRaises(ValueError):
            parse_color_value(123)


class TestLoadCustomPalette(unittest.TestCase):
    """Tests for load_custom_palette in sqanti_browser."""

    def test_load_valid_palette(self):
        from sqanti_browser import load_custom_palette

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(
                {
                    "standard": {"full-splice_match": [107, 174, 214]},
                    "highlight": {},
                },
                f,
            )
            path = f.name
        try:
            std, hl, val = load_custom_palette(path)
            self.assertIn("full-splice_match", std)
            self.assertEqual(std["full-splice_match"], (107, 174, 214))
        finally:
            Path(path).unlink(missing_ok=True)

    def test_load_hex_in_palette(self):
        from sqanti_browser import load_custom_palette

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(
                {"standard": {"full-splice_match": "#6BAED6"}},
                f,
            )
            path = f.name
        try:
            std, _, _ = load_custom_palette(path)
            self.assertEqual(std["full-splice_match"], (107, 174, 214))
        finally:
            Path(path).unlink(missing_ok=True)

    def test_missing_file_raises(self):
        from sqanti_browser import load_custom_palette

        with self.assertRaises(FileNotFoundError):
            load_custom_palette("/nonexistent/palette.json")

    def test_invalid_json_raises(self):
        from sqanti_browser import load_custom_palette

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            f.write("{ invalid json }")
            path = f.name
        try:
            with self.assertRaises(ValueError):
                load_custom_palette(path)
        finally:
            Path(path).unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Filter isoforms
# ---------------------------------------------------------------------------


class TestFilterIsoformsRgbToHex(unittest.TestCase):
    """Tests for filter_isoforms.rgb_to_hex."""

    def test_rgb_to_hex(self):
        from src.filter_isoforms import rgb_to_hex

        self.assertEqual(rgb_to_hex((107, 174, 214)), "#6BAED6")
        self.assertEqual(rgb_to_hex((0, 0, 0)), "#000000")
        self.assertEqual(rgb_to_hex((255, 255, 255)), "#FFFFFF")


class TestFilterIsoformsGetCategoryColor(unittest.TestCase):
    """Tests for filter_isoforms.get_category_color."""

    def test_known_category(self):
        from src.filter_isoforms import get_category_color

        color = get_category_color("full-splice_match")
        self.assertTrue(color.startswith("#") and len(color) == 7)

    def test_unknown_category_default(self):
        from src.filter_isoforms import get_category_color

        color = get_category_color("unknown_cat")
        self.assertEqual(color, "#C8C8C8")  # default gray


class TestFilterIsoformsGenerateCategorySvg(unittest.TestCase):
    """Tests for filter_isoforms.generate_category_svg."""

    def test_returns_svg_string(self):
        from src.filter_isoforms import generate_category_svg

        svg = generate_category_svg("full-splice_match", include_reference=True)
        self.assertIn("<svg", svg)
        self.assertIn("fill=", svg)

    def test_unknown_category_returns_empty(self):
        from src.filter_isoforms import generate_category_svg

        svg = generate_category_svg("unknown_category")
        self.assertEqual(svg, "")


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------


class TestConstants(unittest.TestCase):
    """Tests for src.constants."""

    def test_bed12_base_cols_count(self):
        from src.constants import BED12_BASE_COLS

        self.assertEqual(len(BED12_BASE_COLS), 12)

    def test_ordered_filters_structure(self):
        from src.constants import ORDERED_FILTERS

        for item in ORDERED_FILTERS:
            self.assertIn(item[0], ("text", "range"))
            self.assertGreaterEqual(len(item), 3)

    def test_default_palette_has_known_categories(self):
        from src.constants import DEFAULT_STANDARD_PALETTE

        cats = {"full-splice_match", "novel_in_catalog", "genic_intron"}
        for c in cats:
            self.assertIn(c, DEFAULT_STANDARD_PALETTE)
            self.assertEqual(len(DEFAULT_STANDARD_PALETTE[c]), 3)


# ---------------------------------------------------------------------------
# BedProcessor sort_bed_for_visualization
# ---------------------------------------------------------------------------


class TestBedProcessorSort(unittest.TestCase):
    """Tests for BedProcessor.sort_bed_for_visualization."""

    def _make_mock_converter(self, sort_by="none"):
        class Mock:
            pass

        m = Mock()
        m.sort_by = sort_by
        return m

    def test_sort_by_genomic_position(self):
        from src.bed_processor import BedProcessor

        conv = self._make_mock_converter("none")
        proc = BedProcessor(conv)
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2"],
                "chromStart": [100, 50, 10],
                "associated_transcript": ["A", "A", "B"],
            }
        )
        result = proc.sort_bed_for_visualization(df)
        self.assertEqual(list(result["chrom"]), ["chr1", "chr1", "chr2"])
        self.assertEqual(list(result["chromStart"]), [50, 100, 10])

    def test_sort_with_missing_sort_column(self):
        from src.bed_processor import BedProcessor

        conv = self._make_mock_converter("nonexistent_col")
        proc = BedProcessor(conv)
        df = pd.DataFrame({"chrom": ["chr1"], "chromStart": [100]})
        result = proc.sort_bed_for_visualization(df)
        self.assertEqual(len(result), 1)


# ---------------------------------------------------------------------------
# ValidationTrackBuilder _pack_rgb
# ---------------------------------------------------------------------------


class TestValidationTrackPackRgb(unittest.TestCase):
    """Tests for ValidationTrackBuilder._pack_rgb."""

    def test_pack_rgb(self):
        from src.validation_tracks import ValidationTrackBuilder

        class Mock:
            pass

        b = ValidationTrackBuilder(Mock())
        self.assertEqual(b._pack_rgb(255, 0, 0), 16711680)  # red
        self.assertEqual(b._pack_rgb(0, 0, 0), 0)


if __name__ == "__main__":
    unittest.main()
