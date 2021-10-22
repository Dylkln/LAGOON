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
                        help="The columns names used in attributes")

    return parser.parse_args()


def percentage(part, whole):
  percent = 100 * float(part)/float(whole)
  return str(round(percent, 2)) + "%"


def read_file(file):
    with open(file, "r") as f:
        for line in f:
            yield line


def get_files_from_arg(file):

    files = []
    with open(file, "r") as f_in:
        for line in f_in:
            files.append(line.strip())

    return files


def read_fasta_files(files, output):
    out = open(output, "w")
    for file in files:
        records = SeqIO.parse(open(file), "fasta")
        name = "-".join(file.split("-")[1:2]).split("/")[-1]
        for record in records:
            if "TRINITY" in record.id:
                record.id = record.id.replace("TRINITY", name)
                record.description = record.description.replace("TRINITY", name)
            SeqIO.write(record, out, "fasta")
    out.close()


def create_diamond_db(fasta_file):
    db = "./results/all_data_db"
    try:
        command = ["diamond-aligner", "makedb", "--in", fasta_file, "--db", db]
        subprocess.call(command, stdout=DEVNULL)

    except:
        command = ["diamond", "makedb", "--in", fasta_file, "--db", db]
        subprocess.call(command, stdout=DEVNULL)
    return db


def diamond_blastp(db, fasta_file, output):

    try:
        command = ["diamond-aligner", "blastp", "-d", db, "-q", fasta_file, "-o", output, "-e",
               "1e-5", "--sensitive", "-f", "6", "qseqid", "qlen", "qstart", "qend", "sseqid",
               "slen", "sstart", "send", "length", "pident", "ppos", "score", "evalue", "bitscore"]
        subprocess.call(command, stdout=DEVNULL)

    except:
        command = ["diamond", "blastp", "-d", db, "-q", fasta_file, "-o", output, "-e",
                   "1e-5", "--sensitive", "-f", "6", "qseqid", "qlen", "qstart", "qend", "sseqid",
                   "slen", "sstart", "send", "length", "pident", "ppos", "score", "evalue",
                   "bitscore"]
        subprocess.call(command, stdout=DEVNULL)


def create_vertices_edges_files(ssn):
    command = ["perl", "diamond2graph_V2.pl", "--input", ssn]
    subprocess.call(command, stdout=DEVNULL)


def filter(inputfile, outputfile, cov, ident):
    al_ssn, al_filt, nb_nssn, nb_nfilt = 0, 0, 0, 0
    n_ssn = set([])
    n_filt = set([])

    fieldnames = ["qseqid", "qlen", "qstart", "qend", "sseqid", "slen",
                  "sstart", "send", "length", "pident", "ppos", "score", "evalue",
                  "bitscore"]

    if inputfile.endswith(".gz"):
        i = gzip.open(inputfile, mode="rt")
        f = io.TextIOWrapper(i, encoding="utf-8")
        reader = csv.DictReader(f, delimiter="\t",
                                fieldnames=fieldnames)
    else:
        reader = csv.DictReader(open(inputfile), delimiter="\t",
                                fieldnames=fieldnames)

    writer = csv.DictWriter(open(outputfile, "w"), delimiter="\t",
                            fieldnames=fieldnames)

    for row in reader:

        al_ssn += 1
        n = [row["qseqid"], row["sseqid"]]

        if n[0] not in n_ssn:
            n_ssn.add(n[0])
            nb_nssn += 1

        if not row["pident"] or not row["ppos"]:
            continue

        else:
            if float(row["pident"]) >= ident and float(row["ppos"]) >= cov:
                if n[0] == n[1]:
                    continue
                else:
                    writer.writerow(row)
                    al_filt += 1
                    if n[0] not in n_filt:
                        n_filt.add(n[0])
                        nb_nfilt += 1

    return al_ssn, al_filt, nb_nssn, nb_nfilt


def filter_file(identity, overlap, ssn):
    for i in identity:
        for j in overlap:
            print("*** FILTERING FILE ***")
            print(f"filtering infos ||| coverage : {j}%, identity : {i}%")
            output = f"{ssn}_pcov{int(j)}_pident{int(i)}"
            al_ssn, al_filt, nb_nssn, nb_nfilt = filter(ssn, output, j, i)
            rm = al_ssn - al_filt
            rmnb = nb_nssn - nb_nfilt
            print(f"nb of alignments in base SSN : {al_ssn}")
            print(f"nb of alignments in filtered SSN : {al_filt} ({percentage(al_filt, al_ssn)})")
            print(f"nb of alignments removed : {rm} ({percentage(rm, al_ssn)})")
            print(f"nb of nodes in base SSN : {nb_nssn}")
            print(f"nb of nodes in filtered SSN : {nb_nfilt} ({percentage(nb_nfilt,nb_nssn)})")
            print(f"nb of nodes removed : {rmnb} ({percentage(rmnb, nb_nssn)})")

            output_stats = f"{output}_stats"
            save_stats(output_stats, al_ssn, al_filt, nb_nssn, nb_nfilt)
            print("*** CREATING EDGES AND VERTICES FILE ***")
            diamond2graph(output)


def save_stats(output, al_ssn, al_filt, nb_nssn, nb_nfilt):
    rm = al_ssn - al_filt
    rmnb = nb_nssn - nb_nfilt
    with open(output, "w") as f:
        f.write(f"nb of alignments in base SSN : {al_ssn}\n")
        f.write(f"nb of alignments in filtered SSN : {al_filt} ({percentage(al_filt, al_ssn)})\n")
        f.write(f"nb of alignments removed : {rm} ({percentage(rm, al_ssn)})\n")
        f.write(f"nb of nodes in base SSN : {nb_nssn}\n")
        f.write(f"nb of nodes in filtered SSN : {nb_nfilt} ({percentage(nb_nfilt,nb_nssn)})\n")
        f.write(f"nb of nodes removed : {rmnb} ({percentage(rmnb, nb_nssn)})\n")


def diamond2graph(diamond_output):
    ids = set()
    fieldnames = "from;to;query_length;subject_length;alignment_len;pident;evalue;bitscore"
    edges = open(f"{diamond_output}.edges", "w")
    vertices = open(f"{diamond_output}.vertices", "w")
    edges.write(f"{fieldnames}\n")
    vertices.write("name;prefix\n")
    for line in read_file(diamond_output):
        tab = line.split("\t")
        qsid, qlen, ssid, slen, leng, pid, eval, bscore = tab[0], tab[1], tab[4],\
                                                          tab[5], tab[8], tab[9],\
                                                          tab[12], tab[13]
        edges.write(f"{qsid};{ssid};{qlen};{slen};{leng};{pid};{eval};{bscore}")

        if qsid not in ids:
            ids.add(qsid)
        if ssid not in ids:
            ids.add(ssid)

    ids = sorted(ids)

    for i in ids:
        prefix = i.split(".")[0]
        vertices.write(f"{i};{prefix}\n")


def adapt_row(row, columns):
    d = {}
    n = row[columns[0]]
    c = columns[1:]

    for k, v in row.items():
        if k in c:
            d[k] = v

    for k, v in d.items():
        if v == None:
            d[k] = "NA"

    return d, n, c


def create_attributes_dict(an_files, columns):
    at_dict = {}
    for file in an_files:
        reader = csv.DictReader(open(file), delimiter="\t")
        n = file.split("-")[1]
        for row in reader:
            d, nc, c = adapt_row(row, columns)
            if n not in at_dict.keys():
                at_dict[n] = {}
            if nc not in at_dict[n].keys():
                at_dict[n][nc] = {k : set() for k in d.keys()}
            for k, v in d.items():
                at_dict[n][nc][k].add(v)

    return at_dict


def save_attributes(at_dict, columns):
    path = "./results/attributes/"
    if not os.path.exists(path):
        os.mkdir(path)

    for k, v in at_dict.items():
        tmp = "-".join(k.split("-")[0:2])
        output = f"./results/attributes/{tmp}.attrib"
        fields = ";".join(columns)

        with open(output, "w") as f:
            f.write(f"{fields}\n")
            for k2, v2 in v.items():
                lwrite = [i for i in v2.values()]
                lwrite = ";".join([",".join(i) for i in lwrite])
                f.write(f"{k2};{lwrite}\n")

def main():
    start = time.time()
    args = arguments()


    if args.fasta_files and args.annotation_files and args.columns:
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
#            if not os.path.exists(ssn):
            diamond_blastp(db, fasta, ssn)

        elif args.fasta_files[0].endswith(".pep"):
            fasta = "./results/all_data.fasta"
            read_fasta_files(args.fasta_files, fasta)
            db = create_diamond_db(fasta)
            ssn = "./results/diamond_ssn"
#            if not os.path.exists(ssn):
            diamond_blastp(db, fasta, ssn)

        else:
            files = get_files_from_arg(args.fasta_files[0])
            fasta = "./results/all_data.fasta"
            read_fasta_files(files, fasta)
            db = create_diamond_db(fasta)
            ssn = "./results/diamond_ssn"
#            if not os.path.exists(ssn):
            diamond_blastp(db, fasta, ssn)


    if args.overlap and args.identity:
        print("filtering informations provided")
        ssn = "./results/diamond_ssn"
        filter_file(args.identity, args.overlap, ssn)

    if not args:
        print("Please enter some arguments, see --help / -h for more infos")

    print(f"it took {round((time.time() - start), 2)} seconds")


if __name__ == '__main__':
    main()
