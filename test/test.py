"""This script is used to create all test files used to test all scripts"""


#================================== Modules ==================================#

import argparse

#=============================================================================#

def arguments():
    """
    set arguments
    """

    parser = argparse.ArgumentParser()

    # Mandatory arguments
    parser.add_argument("-f", "--file", dest = "file",
        type = str, required = True,
        help = "the name of the input file")

    return parser.parse_args()

def read_file(file):

    with open(file) as f:
        for line in f:
            yield line.strip()

def get_data(file):

    for index, line in enumerate(read_file(file)):
        llist = line.split()
        if index == 0:
            if llist[0] == llist[1]:
                eq = set(llist)
            elif llist[0] != llist[1]:
                neq = set(llist)
            else:
                eq = set()
                neq = set()
                oth = set(llist)
        else:
            if llist[0] == llist[1]:
                if len(eq) <= 10:
                    eq.add(llist)
            elif llist[0] != llist[1]:
                if len(neq) <= 10:
                    neq.add(llist)
            else:
                if len(oth) <= 50:
                    oth.add(llist)
        if len(eq) == 10 and len(neq) == 10 and len(oth) == 50:
            return eq, neq, oth

    return eq, neq, oth

def main():

    args = arguments()
    eq, neq, oth = get_data(args.file)
