class mcolors:
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def warning(text):
    """
    Write a warning message in bold red
    """
    print(f"{mcolors.RED}{mcolors.BOLD}%s{mcolors.ENDC}" % text)

def message(text):
    """
    Write a normal message
    """
    print(text)

def message_OK(text):
    """
    Write a message in green
    """
    print(f"{mcolors.GREEN}%s{mcolors.ENDC}" % text)

def emphasis(text):
    """
    Write an important message
    """
    print(f"{mcolors.BLUE}{mcolors.BOLD}%s{mcolors.ENDC}" % text)
