import os
import csv
import time

def adapt_row(row, columns, i_dict):

    n = i_dict[row[columns[0]]]
    c = columns[1:]
    d = {k: v for k, v in row.items() if k in c}

    for k, v in d.items():
        if v is None:
            d[k] = "NA"

    return d, n, c


def determine_file(file):
    if type(file) != str:
        file = str(file)
    with open(file) as f:
        line = f.readline()
        print(bool(os.path.exists(line.strip())))
        return bool(os.path.exists(line.strip()))


def get_files_from_arg(file):
    if type(file) != str:
        file = str(file)
    files = []
    with open(file, "r") as f_in:
        for line in f_in:
            files.append(line.strip())

    return files


def treat_annotation(an_files, columns, indices):
    at_dict = {}
    i_dict = {}


    with open(indices, "r") as f:
        for line in f:
            llist = line.split("\t")[:2]
            i_dict[llist[1]] = llist[0]

    for file in an_files:
        reader = csv.DictReader(open(file), delimiter="\t")
        n = file.split("-")[1]
        for row in reader:
            d, nc, c = adapt_row(row, columns, i_dict)
            if n not in at_dict.keys():
                at_dict[n] = {}
            if nc not in at_dict[n].keys():
                at_dict[n][nc] = {k: set() for k in d.keys()}
            for k, v in d.items():
                at_dict[n][nc][k].add(v)
    return at_dict


def create_attributes_dict(an_files, columns, indices):
    if len(an_files) != 1 or not determine_file(an_files):
        return treat_annotation(an_files, columns, indices)
    files = get_files_from_arg(an_files)
    return treat_annotation(files, columns, indices)


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
    with open(str(snakemake.log), "w") as log:
        s = time.time()
        log.write("*** Getting input files and parameters ***\n")
        columns = snakemake.params.columns
        files = snakemake.input.an_files
        indices = snakemake.input.indices
        log.write("*** Getting Attributes from files ***\n")
        at_dict = create_attributes_dict(files, columns, indices)
        log.write("*** Saving Attributes ***\n")
        save_attributes(at_dict, columns)
        e = time.time()
        log.write(f"Operations done in {round(e - s, 2)} seconds")


if __name__ == '__main__':
    main()
