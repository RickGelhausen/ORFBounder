#!/usr/bin/env python3
"""Compatibility CLI for the supported ORFBounder result-table merger.

Sample names come from table columns. Unlike the former helper, this does not
create experiment namespaces from filenames. Use unique sample column names
when combining experiments; repeated names refer to the same sample.
"""

import argparse
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lib.merging import merge_tables


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge ORFBounder Excel or tab-separated CSV/TSV result tables.",
        epilog=("Deprecated entry point: prefer python -m lib.merging -t INPUT... -o OUTPUT.xlsx. "
                "Sample names are read from columns; filenames no longer create experiment namespaces."),
    )
    parser.add_argument("-i", "--input", nargs="+", type=Path, required=True, help="One or more result tables.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output .xlsx path; also writes .csv and .gff.")
    args = parser.parse_args(argv)
    try:
        merge_tables(args.input, args.output)
    except (OSError, ValueError, KeyError) as exc:
        parser.error(f"Cannot merge tables: {exc}")


if __name__ == "__main__":
    main()
