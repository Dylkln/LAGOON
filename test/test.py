"""This script is used to test your program installation and to familiarize yourself with the available
arguments. Please see the list of arguments below."""

# ================================== Modules ================================== #

import argparse
import os
import sys
from Bio import SeqIO
from modules import find as fi
from modules import cat

# ============================================================================= #

def arguments():
    """
    set arguments
    """

    parser = argparse.ArgumentParser(description="This script is used to test your program installation and to \
familiarize yourself with the available arguments. Please see the list of arguments below.")

    # Optional Arguments
    parser.add_argument("-tr", dest="fasta_files",
                        type=str, required=False,
                        help="either a path to find files in fasta format or a file containing the list of all path \
to files in fasta format.")
    parser.add_argument("-an", dest="annotation_files",
                        type=str, required=False,
                        help="either a path to find annotation files or a file containing the list of all path \
to annotation files.")
    parser.add_argument("-ov", "--overlap", dest="overlap",
                        type=float, required=False, default=None, nargs="+",
                        help="the filtration wanted of a diamond output")
    parser.add_argument("-id", "--identity", dest="identity",
                        type=float, required=False, default=None, nargs="+",
                        help="the filtration wanted of a diamond output")

    return parser.parse_args()


def read_fasta_files(files, output):
    for file in files
        f_seq = SeqIO.parse(open(file), "fasta")
        with open(output, "w") as out:
            for s in f_seq:
                out.write(s)



def read_file(file):
    with open(file, "r") as f:
        for line in f:
            yield line.strip()


def main():
    args = arguments()
    if not args:
        print("Please enter some arguments, see --help / -h for more infos")

    if args.path and args.pattern:
        files = fi.find_file(args.pattern, args.path)
        path = "./results/"
        if not os.path.exists(path):
            os.mkdir(path)
    elif args.overlap and args.identity:
        pass
    else:
        pass


if __name__ == '__main__':
    main()
