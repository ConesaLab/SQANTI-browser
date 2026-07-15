"""Install the bundled SQANTI-browser Agent Skill into a Claude skills directory.

Exposed as the console script ``sqanti-browser-install-skills`` (see pyproject.toml).

Usage
-----
    sqanti-browser-install-skills            # copy into ~/.claude/skills/sqanti-browser
    sqanti-browser-install-skills --force    # overwrite an existing install
    sqanti-browser-install-skills --print-path   # print the bundled skill dir, don't install

The last form is handy for pointing a custom skills location at the bundled files
without copying, e.g.::

    export CLAUDE_SKILLS_PATH="$(sqanti-browser-install-skills --print-path)"
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

SKILL_NAME = "sqanti-browser"
_IGNORE = shutil.ignore_patterns("__pycache__", "*.pyc", ".ipynb_checkpoints", ".DS_Store")


def _bundled_skill_dir() -> Path:
    """Directory that holds the skill files (SKILL.md + references + assets)."""
    return Path(__file__).resolve().parent / "data"


def _default_dest() -> Path:
    """Default install location for Claude Code user skills."""
    return Path.home() / ".claude" / "skills" / SKILL_NAME


def install_skill(dest: Path | None = None, force: bool = False) -> Path:
    """Copy the bundled skill into ``dest`` (default ~/.claude/skills/sqanti-browser).

    Returns the destination path. Raises on a missing source or an existing
    destination when ``force`` is False.
    """
    src = _bundled_skill_dir()
    if not (src / "SKILL.md").is_file():
        raise FileNotFoundError(
            f"Bundled skill is missing SKILL.md (looked in {src}). "
            "Is the package installed with its data files?"
        )

    dest = dest or _default_dest()
    if dest.exists():
        if not force:
            raise FileExistsError(
                f"{dest} already exists. Re-run with --force to overwrite "
                "(e.g. after upgrading SQANTI-browser)."
            )
        shutil.rmtree(dest)

    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dest, ignore=_IGNORE)
    return dest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="sqanti-browser-install-skills",
        description="Install the SQANTI-browser Claude Agent Skill into your skills directory.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing installation (use after upgrading SQANTI-browser).",
    )
    parser.add_argument(
        "--print-path",
        action="store_true",
        help="Print the bundled skill directory and exit without installing.",
    )
    parser.add_argument(
        "--dest",
        type=Path,
        default=None,
        help=f"Install destination (default: {_default_dest()}).",
    )
    args = parser.parse_args(argv)

    if args.print_path:
        print(_bundled_skill_dir())
        return 0

    try:
        dest = install_skill(dest=args.dest, force=args.force)
    except (FileNotFoundError, FileExistsError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"Installed the SQANTI-browser skill to {dest}")
    print("Restart your agent (or start a new session) to pick it up.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
