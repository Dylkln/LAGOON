import argparse
import csv
import io
import gzip


def arguments():
    """
    set arguments
    """

    parser = argparse.ArgumentParser()

    # Optional Arguments
    parser.add_argument("-f", "--file", dest="file",
                        type=str, required=False, default=None,
                        help="the ssn file")
    parser.add_argument("-ov", "--overlap", dest="overlap",
                        type=float, required=False, default=None, nargs="+",
                        help="the filtration wanted of a diamond output")
    parser.add_argument("-id", "--identity", dest="identity",
                        type=float, required=False, default=None, nargs="+",
                        help="the filtration wanted of a diamond output")
    return parser.parse_args()


def percentage(part, whole):
    percent = 100 * float(part) / float(whole)
    return str(round(percent, 2)) + "%"


def filter(inputfile, outputfile, cov, ident):
    al_ssn, al_filt, nb_nssn, nb_nfilt = 0, 0, 0, 0
    n_ssn = set()
    n_filt = set()

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

        if (
                row["pident"]
                and row["ppos"]
                and float(row["pident"]) >= ident
                and float(row["ppos"]) >= cov
                and n[0] != n[1]
        ):
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
            print(f"nb of nodes in filtered SSN : {nb_nfilt} ({percentage(nb_nfilt, nb_nssn)})")
            print(f"nb of nodes removed : {rmnb} ({percentage(rmnb, nb_nssn)})")

            output_stats = f"{output}_stats"
            save_stats(output_stats, al_ssn, al_filt, nb_nssn, nb_nfilt)


def save_stats(output, al_ssn, al_filt, nb_nssn, nb_nfilt):
    rm = al_ssn - al_filt
    rmnb = nb_nssn - nb_nfilt
    with open(output, "w") as f:
        f.write(f"nb of alignments in base SSN : {al_ssn}\n")
        f.write(f"nb of alignments in filtered SSN : {al_filt} ({percentage(al_filt, al_ssn)})\n")
        f.write(f"nb of alignments removed : {rm} ({percentage(rm, al_ssn)})\n")
        f.write(f"nb of nodes in base SSN : {nb_nssn}\n")
        f.write(f"nb of nodes in filtered SSN : {nb_nfilt} ({percentage(nb_nfilt, nb_nssn)})\n")
        f.write(f"nb of nodes removed : {rmnb} ({percentage(rmnb, nb_nssn)})\n")


# def diamond2graph(diamond_output):
#    ids = set()
#    fieldnames = "from;to;query_length;subject_length;alignment_len;pident;evalue;bitscore"
#    edges = open(f"{diamond_output}.edges", "w")
#    vertices = open(f"{diamond_output}.vertices", "w")
#    edges.write(f"{fieldnames}\n")
#    vertices.write("name;prefix\n")
#    for line in read_file(diamond_output):
#        tab = line.split("\t")
#        qsid, qlen, ssid, slen, leng, pid, eval, bscore = tab[0], tab[1], tab[4],\
#                                                          tab[5], tab[8], tab[9],\
#                                                          tab[12], tab[13]
#        edges.write(f"{qsid};{ssid};{qlen};{slen};{leng};{pid};{eval};{bscore}")

#        if qsid not in ids:
#            ids.add(qsid)
#        if ssid not in ids:
#            ids.add(ssid)

#   ids = sorted(ids)

#    for i in ids:
#        prefix = i.split(".")[0]
#        vertices.write(f"{i};{prefix}\n")


def main():
    args = arguments()
    print("filtering informations provided")
    filter_file(args.identity,
                args.overlap,
                args.file)


if __name__ == '__main__':
    main()
