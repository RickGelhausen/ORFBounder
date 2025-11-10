"""
Messaging functions for terminal output with colors.
"""
import sys

# ANSI color codes
BLUE = '\033[94m'
CYAN = '\033[96m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'


def error(text):
    """Write an error message in bold red to stderr."""
    print(f"{RED}{BOLD}{text}{ENDC}", file=sys.stderr)


def error_list(prefix_text, suffix_text, input_description, expected_list, input_list):
    """
    Display a colored list showing missing (red) and present (green) items.

    Args:
        prefix_text: Text before the list
        suffix_text: Text after the list
        input_description: Description of the input
        expected_list: List of expected items
        input_list: List of actual items
    """
    missing_entries = set(expected_list) - set(input_list)
    indent = " " * (len(input_description) + 1)

    item_list = []
    for i, item in enumerate(expected_list):
        color = RED if item in missing_entries else GREEN
        prefix = " " if i == 0 else indent
        item_list.append(f"{prefix}{color}{BOLD}{item}{ENDC}\n")

    error_msg = (
        f"{RED}{BOLD}{prefix_text}{ENDC}"
        f"{RED}{BOLD}{input_description}{ENDC}"
        f"{''.join(item_list)}"
        f"{RED}{BOLD}{suffix_text}{ENDC}"
    )

    print(error_msg, file=sys.stderr)


def success(text):
    """Write a success message in green."""
    print(f"{GREEN}{text}{ENDC}")


def message(text):
    """Write a normal message in blue."""
    print(f"{BLUE}{text}{ENDC}")


def warning(text):
    """Write a warning message in yellow to stderr."""
    print(f"{YELLOW}{text}{ENDC}", file=sys.stderr)