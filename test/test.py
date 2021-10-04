"""This script is used to create all test files used to test all scripts"""

# ================================== Modules ==================================#

import argparse


# =============================================================================#

def arguments():
    """
    set arguments
    """

    parser = argparse.ArgumentParser(description="This script is used to test \
your program installation and to familiarize yourself with the available \
arguments. Please see the list of arguments below.")

    # Optional Arguments
    parser.add_argument("-pt", "--pattern", dest="pattern",
                        type=str, required=False,
                        help="the pattern of the files to be found")
    parser.add_argument("-ph", "--path", dest="path",
                        type=str, required=False,
                        help="absolute path to directory for the pattern")

    return parser.parse_args()


def read_file(file):
    with open(file) as f:
        for line in f:
            yield line.strip()


def main():
    args = arguments()
    if not args:
        print("Please enter some arguments, see --help for more infos")


if __name__ == '__main__':
    main()
