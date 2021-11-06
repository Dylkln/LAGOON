import os
import csv
import argparse


def arguments():
    """
    set arguments
    """

    parser = argparse.ArgumentParser()

    # Optional Arguments
    parser.add_argument("-an", dest="annotation_files",
                        type=str, required=False, default=None, nargs="+",
                        help="either a path to find annotation files or a file containing the "
                             "list of all path to annotation files.")
    parser.add_argument("-cn", "--columns", dest="columns",
                        type=str, required=False, default=None, nargs="+",
                        help="the columns names used in attributes")
    return parser.parse_args()


def adapt_row(row, columns):
    n = row[columns[0]]
    c = columns[1:]
    d = {k: v for k, v in row.items() if k in c}

    for k, v in d.items():
        if v is None:
            d[k] = "NA"

    return d, n, c


def determine_file(file):
    with open(file) as f:
        line = f.readline()
        return bool(os.path.exists(line.strip()))


def get_files_from_arg(file):
    files = []
    with open(file, "r") as f_in:
        for line in f_in:
            files.append(line.strip())

    return files


def treat_annotation(an_files, columns):
    at_dict = {}
    for file in an_files:
        reader = csv.DictReader(open(file), delimiter="\t")
        n = file.split("-")[1]
        for row in reader:
            d, nc, c = adapt_row(row, columns)
            if n not in at_dict.keys():
                at_dict[n] = {}
            if nc not in at_dict[n].keys():
                at_dict[n][nc] = {k: set() for k in d.keys()}
            for k, v in d.items():
                at_dict[n][nc][k].add(v)
    return at_dict


def create_attributes_dict(an_files, columns):
    if len(an_files) != 1 or determine_file(an_files):
        return treat_annotation(an_files, columns)
    files = get_files_from_arg(an_files)
    return treat_annotation(files, columns)


def save_attributes(at_dict, columns):
    path = "./results/attributes/"
    if not os.path.exists(path):
        os.mkdir(path)

    for k, v in at_dict.items():
        output = f"./results/attributes/{k}.attributes"
        with open(output, "w") as f:
            f.write(";".join(columns) + "\n")
            for k2, v2 in v.items():
                lwrite = [i for i in v2.values()]
                lwrite = ";".join(",".join(i) for i in lwrite)
                f.write(f"{k2};{lwrite}\n")


def main():
    args = arguments()
    at_dict = create_attributes_dict(args.annotation_files, args.columns)
    save_attributes(at_dict, args.columns)


if __name__ == '__main__':
    main()
