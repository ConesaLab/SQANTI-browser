# Changelog

All notable changes to SQANTI-browser are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2025-01-19

### Added

- **`--category-tracks`**: Create tracks only for specified structural categories. Comma-separated abbreviated names (FSM, ISM, NIC, NNC, antisense, genic_intron, genic_genomic, intergenic, fusion). Example: `--category-tracks FSM,ISM,NIC`
- **Packaging**: `pyproject.toml` for modern Python packaging
- **CLI entry point**: `sqanti_browser` command (after `pip install -e .`)
- **Modular architecture**: Split core logic into `BedProcessor`, `HubGenerator`, `ValidationTrackBuilder`
- **Type hints**: Added to public APIs across the codebase
- **Unit tests**: `tests/test_unit.py` for utilities, color parsing, constants
- **Edge case tests**: `tests/test_edge_cases.py` for missing files, empty inputs, CLI handling
- **Constants module**: Centralized filter limits (e.g. exon count ≤400, transcript length ≤150 kb) in `src/constants.py`
- **Documentation**: "How to run" clarification (`python -m sqanti_browser` vs `sqanti_browser` vs `python sqanti_browser.py`)

### Changed

- **Renamed**: `sqanti3_to_UCSC.py` → `sqanti_browser.py` (git history preserved)
- **Recommended run command**: Use `python -m sqanti_browser` so the active Python (e.g. conda) is always used
- **Project layout**: `src/`, `tests/`, `example/` directories; `filter_isoforms.py` moved to `src/`
- **CI**: GitHub Actions updated with `pip install -e .`, `PYTHONPATH`, unit and edge tests
- **Filter limits**: Increased max exon count (200→400) and transcript length (100 kb→150 kb) for biologically realistic genes (TTN, DMD)

### Fixed

- CI `cwd` and `PYTHONPATH` for reliable test execution
- Import paths after module reorganization

## [1.0.0] - (prior releases)

Initial functionality:

- Convert SQANTI3 GTF and classification to UCSC track hubs
- Category-specific tracks, validation tracks (CAGE, PolyA, STAR, reference)
- Trix search, interactive HTML tables, custom palettes
- SQANTI-reads integration
