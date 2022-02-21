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

def error_list(prefix_text, suffix_text, input_description, expected_list, input_list):
    """
    Write an error message in bold red,
    """
    missing_entries = list(set(expected_list) - set(input_list))
    description_length = len(input_description)

    item_list = []
    for item in expected_list:
        if len(item_list) == 0:
            if item in missing_entries:
                item_list.append(" "+f"{mcolors.RED}{mcolors.BOLD}%s{mcolors.ENDC}\n" % item)
            else:
                item_list.append(" "+f"{mcolors.GREEN}{mcolors.BOLD}%s{mcolors.ENDC}\n" % item)
        else:
            if item in missing_entries:
                item_list.append(" "*(description_length+1)+ f"{mcolors.RED}{mcolors.BOLD}%s{mcolors.ENDC}\n" % item)
            else:
                item_list.append(" "*(description_length+1)+ f"{mcolors.GREEN}{mcolors.BOLD}%s{mcolors.ENDC}\n" % item)

    print(f"{mcolors.RED}{mcolors.BOLD}%s{mcolors.ENDC}" % prefix_text\
         +f"{mcolors.RED}{mcolors.BOLD}%s{mcolors.ENDC}" % input_description\
         +"".join(item_list)\
         +f"{mcolors.RED}{mcolors.BOLD}%s{mcolors.ENDC}" % suffix_text)

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
