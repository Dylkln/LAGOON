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


# ============================================================================= #

def arguments():
    """
    set arguments
    """

    parser = argparse.ArgumentParser(description="This script is used to test your program installation and to \
familiarize yourself with the available arguments. Please see the list of arguments below.")

    # Optional Arguments
    parser.add_argument("-tr", dest="fasta_files",
                        type=str, required=False, nargs="+",
                        help="either a path to find files in fasta format or a file containing the list of all path \
to files in fasta format.")
    parser.add_argument("-an", dest="annotation_files",
                        type=str, required=False, nargs="+",
                        help="either a path to find annotation files or a file containing the list of all path \
to annotation files.")
    parser.add_argument("-ov", "--overlap", dest="overlap",
                        type=float, required=False, default=None, nargs="+",
                        help="the filtration wanted of a diamond output")
    parser.add_argument("-id", "--identity", dest="identity",
                        type=float, required=False, default=None, nargs="+",
                        help="the filtration wanted of a diamond output")

    return parser.parse_args()


def percentage(part, whole):
  percent = 100 * float(part)/float(whole)
  return str(round(percent, 2)) + "%"


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
    command = ["diamond-aligner", "makedb", "--in", fasta_file, "--db", db]
    subprocess.call(command, stdout=DEVNULL)
    return db


def diamond_blastp(db, fasta_file, output):
    command = ["diamond-aligner", "blastp", "-d", db, "-q", fasta_file, "-o", output, "-e", "1e-5", "--sensitive", "-f",
               "6", "qseqid", "qlen", "qstart", "qend", "sseqid", "slen", "sstart", "send", "length", "pident", "ppos",
               "score", "evalue", "bitscore"]
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
            create_vertices_edges_files(output)


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


def main():
    args = arguments()


    if args.fasta_files and args.annotation_files:
        path = "./results/"
        if not os.path.exists(path):
            os.mkdir(path)
        fasta = "./results/all_data.fasta"
        read_fasta_files(args.fasta_files, fasta)
        db = create_diamond_db(fasta)
        ssn = "./results/diamond_ssn"
        if not os.path.exists(ssn):
            diamond_blastp(db, fasta, ssn)


    if args.overlap and args.identity:
        print("filtering informations provided")
        ssn = "./results/diamond_ssn"
        filter_file(args.identity, args.overlap, ssn)


    else:
        print("Please enter some arguments, see --help / -h for more infos")


if __name__ == '__main__':
    main()
