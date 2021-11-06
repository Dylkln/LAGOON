from Bio import SeqIO


def read_fasta_files(files, output):
    with open(output, "w") as out:
        for file in files:
            records = SeqIO.parse(open(file), "fasta")
            name = "-".join(file.split("-")[1:2]).split("/")[-1]
            for record in records:
                if "TRINITY" in record.id:
                    record.id = record.id.replace("TRINITY", name)
                    record.description = record.description.replace("TRINITY", name)
                SeqIO.write(record, out, "fasta")


def main():
    files = snakemake.input
    output = str(snakemake.output)
    read_fasta_files(files, output)


if __name__ == '__main__':
    main()
