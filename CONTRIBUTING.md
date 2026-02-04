# Contributing to SQANTI-browser

Thank you for your interest in contributing! This document explains how to set up a development environment and submit changes.

## Development Setup

### 1. Clone the repository

```bash
git clone https://github.com/ConesaLab/SQANTI-browser.git
cd SQANTI-browser
```

### 2. Install UCSC tools

```bash
bash install_ucsc_tools.sh
```

Or with conda:

```bash
conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed ucsc-bedtobigbed ucsc-ixixx ucsc-hubcheck ucsc-twobitinfo ucsc-fatotwobit
```

### 3. Install the package in editable mode

```bash
pip install -e .
```

### 4. (Optional) Install dev dependencies

```bash
pip install -e ".[dev]"
```

## Running Tests

From the project root:

```bash
# Integration tests (requires UCSC tools and example data)
python tests/test_sqanti_browser.py

# Quick environment check (no conversion)
python tests/test_sqanti_browser.py --install-only

# Unit tests (no UCSC tools needed)
python tests/test_unit.py -v

# Edge case tests
python tests/test_edge_cases.py -v
```

All tests should pass before submitting a pull request.

## Code Style

- Use type hints on public functions
- Follow existing formatting (PEP 8)
- Prefer `pathlib.Path` over string paths
- Add docstrings for new public APIs

## Submitting Changes

1. **Fork** the repository on GitHub
2. **Create a branch** from `main`: `git checkout -b feature/your-feature-name` or `fix/your-bug-fix`
3. **Make your changes**, and ensure tests pass
4. **Commit** with a clear message: `git commit -m "Add feature X"` or `git commit -m "Fix Y"`
5. **Push** to your fork: `git push origin feature/your-feature-name`
6. **Open a pull request** against `main`

### Pull Request Guidelines

- Describe what the change does and why
- Reference any related issues (e.g. "Fixes #123")
- Ensure CI passes (GitHub Actions runs automatically)

## Documentation

- User docs live in `wiki/`. Edit markdown files there; they are published to the GitHub wiki.
- Update `CHANGELOG.md` for user-facing changes under `[Unreleased]`.

## Questions?

- Open an [issue](https://github.com/ConesaLab/SQANTI-browser/issues) for questions or bugs
- See the [wiki](https://github.com/ConesaLab/SQANTI-browser/wiki) for usage documentation
