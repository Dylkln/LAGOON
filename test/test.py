"""This script is used to test your program installation and to familiarize yourself with the available
arguments. Please see the list of arguments below."""

# ================================== Modules ================================== #


import argparse
import os
import sys
import subprocess
from subprocess import DEVNULL
from Bio import SeqIO
import csv
import io
import gzip
import glob
import time


# ============================================================================= #

def arguments():
    """
    set arguments
    """

    parser = argparse.ArgumentParser(description="This script is used to test your program installation and to \
familiarize yourself with the available arguments. Please see the list of arguments below.")

    # Optional Arguments
    parser.add_argument("-tr", dest="fasta_files",
                        type=str, required=False, default=None, nargs="+",
                        help="either a path to find files in fasta format or a file containing the list of all path \
to files in fasta format.")
    parser.add_argument("-an", dest="annotation_files",
                        type=str, required=False, default=None, nargs="+",
                        help="either a path to find annotation files or a file containing the list of all path \
to annotation files.")
    parser.add_argument("-ov", "--overlap", dest="overlap",
                        type=float, required=False, default=None, nargs="+",
                        help="the filtration wanted of a diamond output")
    parser.add_argument("-id", "--identity", dest="identity",
                        type=float, required=False, default=None, nargs="+",
                        help="the filtration wanted of a diamond output")
    parser.add_argument("-cn", "--columns", dest="columns",
                        type=str, required=False, default=None, nargs="+",
                        help="the columns names used in attributes")
    parser.add_argument("-g", "--graph", dest="graph",
                        type=str, required=False, default=None, nargs="+",
                        help="a sequence similarity network graph")
    parser.add_argument("-t", "--table", dest="table",
                        type=str, required=False, default=None, nargs="+",
                        help="a table used to create the sequence similarity network graph")

    return parser.parse_args()

def get_index(files):

    nodes = set()
    n_dict = {}
    i = 0

    for file in files:
        for line in read_file(file):
            l = line.split("\t")[0]
            if l not in nodes:
                nodes.add(l)
                n_dict[l] = i
                i += 1

        for line in read_file(file):
            l = line.split("\t")[4]
            if l not in nodes:
                nodes.add(l)
                n_dict[l] = i
                i += 1

    return n_dict


def main():
    start = time.time()
    args = arguments()
    x = None

    if args.fasta_files and args.annotation_files and args.columns:
        x = 1
        path = "./results/"
        if not os.path.exists(path):
            os.mkdir(path)


        at_dict = create_attributes_dict(args.annotation_files, args.columns)
        save_attributes(at_dict, args.columns)

        if len(args.fasta_files) > 1:
            fasta = "./results/all_data.fasta"
            read_fasta_files(args.fasta_files, fasta)
            db = create_diamond_db(fasta)
            ssn = "./results/diamond_ssn"
            if not os.path.exists(ssn):
                diamond_blastp(db, fasta, ssn)

        elif args.fasta_files[0].endswith(".pep"):
            fasta = "./results/all_data.fasta"
            read_fasta_files(args.fasta_files, fasta)
            db = create_diamond_db(fasta)
            ssn = "./results/diamond_ssn"
            if not os.path.exists(ssn):
                diamond_blastp(db, fasta, ssn)

        else:
            files = get_files_from_arg(args.fasta_files[0])
            fasta = "./results/all_data.fasta"
            read_fasta_files(files, fasta)
            db = create_diamond_db(fasta)
            ssn = "./results/diamond_ssn"
            if not os.path.exists(ssn):
                diamond_blastp(db, fasta, ssn)


    if args.overlap and args.identity:
        x = 1
        print("filtering informations provided")
        ssn = "./results/diamond_ssn"
        filter_file(args.identity, args.overlap, ssn)

    if not x:
        print("Please enter some arguments, see --help / -h for more infos")

    print(f"it took {round((time.time() - start), 2)} seconds")


if __name__ == '__main__':
    main()
