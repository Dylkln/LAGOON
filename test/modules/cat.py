import os
import time
from Bio import SeqIO


def detect_files(files):

    fs = []
    for file in files:
        with open(file, "r") as f:
            if os.path.exists(f.readline()):
                for line in f:
                    fs.append(line.strip())

    if fs:
        return fs
    return files


def read_fasta_files(files, output):

    i = 1

    with open(output[0], "w") as out:
        with open(output[1], "w") as index:
            for file in files:
                records = SeqIO.parse(open(file), "fasta")
                name = "-".join(file.split("-")[1:2]).split("/")[-1]
                name = name.replace("_", "-") + "-1"
                for record in records:
                    if "TRINITY" in record.id:
                        record.id = record.id.replace("TRINITY", name)
                        record.description = record.description.replace("TRINITY", name)
                    record.id = record.id.replace("_", "-")
                    record.description = record.description.replace("_", "-")
                    index.write(f"{i}\t{record.id}\t{record.description}\n")
                    record.id = str(i)
                    record.description = ""
                    SeqIO.write(record, out, "fasta")
                    i += 1


def main():
    with open(str(snakemake.log), "w") as log:
        s = time.time()
        log.write("*** Getting input and output files ***\n")
        files = detect_files(snakemake.input)
        output = list(snakemake.output)
        log.write("*** Concatenation ***\n")
        read_fasta_files(files, output)
        e = time.time()
        log.write(f"Concatenation of your files done in {round(e - s, 2)} seconds")


if __name__ == '__main__':
    main()
