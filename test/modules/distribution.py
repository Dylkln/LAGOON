import csv
import io
import gzip
import time
import re
import matplotlib.pyplot as plt
from collections import Counter


def read_diamond_output(file):
    file = str(file)
    fieldnames = ["qseqid", "qlen", "qstart", "qend", "sseqid", "slen", "sstart",
                  "send", "length", "pident", "ppos", "score", "evalue",
                  "bitscore"]

    if file.endswith(".gz"):
        i = gzip.open(file, mode="rt")
        f = io.TextIOWrapper(i, encoding="utf-8")
        return csv.DictReader(f, delimiter="\t",
                              fieldnames=fieldnames)
    return csv.DictReader(open(file), delimiter="\t",
                          fieldnames=fieldnames)


def fill_value_dict(reader):
    values_dict = {"pident": [], "ppos": [], "evalue": []}
    for row in reader:
        for k, v in row.items():
            if k in values_dict and k != "evalue":
                values_dict[k].append(float(v))
            if k in values_dict and k == "evalue":
                values_dict[k].append(v)
    return values_dict


def find_graph_file(outputs, key):
    d = {"pident" : "identity", "ppos" : "overlap", "evalue" : "evalue"}
    for k, v in d.items():
        if k == key:
            for file in outputs:
                if re.search(v, file):
                    return file


def save_list_in_file(l, k):
    with open(f"{k}_values", "w") as f:
        for i in l:
            print(i, file=f)


def main():
    with open(str(snakemake.log), "w") as log:
        s = time.time()
        log.write("*** Getting values from alignments ***\n")
        reader = read_diamond_output(snakemake.input)
        values_dict = fill_value_dict(reader)
        for k, v in values_dict.items():
            tmp = find_graph_file(snakemake.output, k)
            g = plt.figure()
            if k != "evalue":
                c = Counter(v)
            else:
                eval = [int(i.split("e")[1]) for i in v]
                c = Counter(eval)
            plt.bar(c.keys(), c.values())
            plt.savefig(tmp)
        e = time.time()
        log.write(f"Operations done in {round(e - s, 2)} seconds")

if __name__ == '__main__':
    main()