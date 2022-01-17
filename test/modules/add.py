# ================================== Modules ==================================#

import re
import csv
import time

# =============================================================================#


def read_file(file):
    """
    Reads a CSV file with a comma delimiter

    Parameters
    ----------

        file : a CSV file

    Yields
    -------

        row : a row of the CSV file
        col : the fieldnames of the file
    """

    with open(file, "r") as f_in:
        col = f_in.readline()
        next(f_in)
        for line in f_in:
            yield col, line.strip()


def read_csv(file):
    """
    Reads a file as CSV

    Parameters
    ----------

        file : a CSV file

    Yields
    -------

        reader : a CSV.Dictreader object
    """
    with open(file, "r") as f_in:
        dialect = csv.Sniffer().sniff(f_in.readline())
        f_in.seek(0)
        reader = csv.DictReader(f_in, dialect=dialect)
        yield from reader


def get_files_from_argument(file):
    """
    Retrieves a list of filenames

    Parameters
    ----------

        file : a file given with the argument -f

    Returns
    -------

        files : a list of files contained in the file
    """

    files = []
    with open(file, "r") as f_in:
        for line in f_in:
            files.append(line.strip())

    return files


def get_fname(names, files, i_dict):
    """
    Retrieves the file name containing the attributes of an ORF ID

    Parameters
    ----------

        n : an ORF ID
        files : a list of filenames

    Returns
    -------

        file : the file containing the attributes of the ORF ID
    """
    fn_dict = {}
    for n in names:
        iname = i_dict[str(n)]
        for file in files:
            if file not in fn_dict:
                fn_dict[file] = set()
            i = "_".join(iname.split("-")[:2])
            if re.search(i, file):
                fn_dict[file].add(n)

    return fn_dict


def get_prefix(n):
    """
    Retrieves the ORF prefix based on an ORF ID

    Parameters
    ----------

        n : an ORF ID

    Returns
    -------

        an ORF ID prefix
    """

    pr = n.replace("-", ".").split(".")[:5]

    if pr[3] != "Transcript":
        del pr[-1]
    return "-".join(pr)


def get_rows(indices, file, columns):
    """
    Retrieves a list of rows containing the ORF name, prefix and attributes

    Parameters
    ----------

        columns : columns names
        nset : a set of ORF IDs
        file : the file where to search the ORFs

    Returns
    -------

        rlist : a list of rows
    """

    rlist = []

    for row in read_csv(file):
        if int(row[columns[0]]) in indices:
            row["name"] = row[columns[0]]
            row.pop(columns[0])
            rlist.append(row)

    return rlist


def write_rows(writer, rows):
    """
    Writes all rows contained in a list of rows in a file

    Parameters
    ----------

        writer : a csv.DictWriter
        rows : a list of rows
    """

    for row in rows:
        writer.writerow(row)


def find_output(file, outputs):
    """
    Retrieves the name of the file needed in a list of files

    Parameters
    ----------

        file : a file to search for
        outputs: a list of files

    Returns
    -------

        f : a file name
    """
    for f in outputs:
        name = file.split(".")[0]
        if name in f:
            return f


def main():
    # Opens Log file
    with open(str(snakemake.log), "w") as log:
        s = time.time()
        log.write("*** Getting input files and parameters ***\n")

        # Get fieldnames
        fieldnames = ["name"] + snakemake.params.columns[1:]

        # get index
        i_dict = {}
        with open(snakemake.input.indices, "r") as f:
            for line in f:
                llist = line.split("\t")[:2]
                i_dict[llist[0]] = llist[1]

        log.write("*** Adding Attributes to the output file ***\n")

        # add attributes to file
        for file in snakemake.input.vertices:
            out = find_output(file, snakemake.output)
            f_out = open(out, "w")
            writer = csv.DictWriter(f_out, delimiter=";", fieldnames=fieldnames)
            writer.writeheader()
            name_set = set()

            for col, line in read_file(file):
                if line not in name_set:
                    name_set.add(int(line))

            fn_dict = get_fname(name_set, snakemake.input.attrib, i_dict)

            for k, v in fn_dict.items():
                rows = get_rows(v, k, snakemake.params.columns)
                write_rows(writer, rows)

        e = time.time()
        log.write(f"Operations done in {round(e - s, 2)} seconds")


if __name__ == '__main__':
    main()
