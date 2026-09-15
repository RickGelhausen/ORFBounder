"""
Messaging functions for terminal output with colors.
"""
import os
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


def _color_enabled(stream) -> bool:
    """Use ANSI styling only for an interactive stream that permits color."""
    if "NO_COLOR" in os.environ or os.environ.get("TERM") == "dumb":
        return False
    try:
        return bool(stream.isatty())
    except (AttributeError, OSError):
        return False


def _styled(text, *styles, stream):
    text = str(text)
    return f"{''.join(styles)}{text}{ENDC}" if _color_enabled(stream) else text


def error(text):
    """Write an error message in bold red to stderr."""
    print(_styled(text, RED, BOLD, stream=sys.stderr), file=sys.stderr)


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
        item_list.append(f"{prefix}{_styled(item, color, BOLD, stream=sys.stderr)}\n")

    error_msg = (
        f"{_styled(prefix_text, RED, BOLD, stream=sys.stderr)}"
        f"{_styled(input_description, RED, BOLD, stream=sys.stderr)}"
        f"{''.join(item_list)}"
        f"{_styled(suffix_text, RED, BOLD, stream=sys.stderr)}"
    )

    print(error_msg, file=sys.stderr)


def success(text):
    """Write a success message in green."""
    print(_styled(text, GREEN, stream=sys.stdout))


def message(text):
    """Write a normal message in blue."""
    print(_styled(text, BLUE, stream=sys.stdout))


def warning(text):
    """Write a warning message in yellow to stderr."""
    print(_styled(text, YELLOW, stream=sys.stderr), file=sys.stderr)
