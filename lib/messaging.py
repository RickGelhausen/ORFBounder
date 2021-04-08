import sys

class mcolors:
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def error(text):
    """
    Write an error message in bold red, leading to a crash and terminate the program.
    """
    print(f"{mcolors.RED}{mcolors.BOLD}%s{mcolors.ENDC}" % text)
    sys.exit()

def success(text):
    """
    Write a message in green
    """
    print(f"{mcolors.GREEN}%s{mcolors.ENDC}" % text)

def message(text):
    """
    Write a normal message
    """
    print(f"{mcolors.BLUE}%s{mcolors.ENDC}" % text)

def warning(text):
    """
    Write a warning message in yellow.
    """
    print(f"{mcolors.YELLOW}%s{mcolors.ENDC}" % text)
