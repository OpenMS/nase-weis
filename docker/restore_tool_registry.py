"""Restore an omitted registry from the matching OpenMS tool declarations.

The current upstream full image copies the source share folder but misses
build-generated TOOLS/*.tsv. Only parse names/categories, never execute CMake.
"""

import argparse
import re
from pathlib import Path


def restore_registry(definitions: Path, tool_dir: Path, output: Path) -> int:
    """Register the declared TOPP tools actually present in this runtime."""
    declarations = re.findall(
        r'openms_topp_tool\(\s*(\w+)\s+"([^"\n]+)"',
        definitions.read_text(encoding="utf-8"),
    )
    entries = sorted(
        {
            (name, category)
            for name, category in declarations
            if (tool_dir / name).is_file()
        }
    )
    if not {"DecoyDatabase", "NucleicAcidSearchEngine"}.issubset(
        {name for name, _ in entries}
    ):
        raise ValueError(
            "The supplied OpenMS declarations/runtime do not contain both NASEWEIS tools."
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        "# Generated from the matching OpenMS release declarations.\n"
        + "".join(f"{name}\t{category}\n" for name, category in entries),
        encoding="utf-8",
    )
    return len(entries)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("definitions", type=Path)
    parser.add_argument("tool_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(
        f"Restored {restore_registry(args.definitions, args.tool_dir, args.output)} TOPP registrations."
    )
