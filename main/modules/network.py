# =========================================== Modules ============================================ #

from collections import Counter
import igraph as ig
import pandas as pd
import csv


# ================================================================================================ #

def percentage(part, whole):
    """
    Retrieves the percentage
    """
    percent = 100 * float(part) / float(whole)
    return round(percent, 2)


def db(string):
    """
    Retrieves the name of the database based on the string provided in argument
    """
    lst_of_lst = [["PFAM", "PF"], ["SMART", "SM"], ["PROSITEPROFILES", "PS5"],
                  ["GENE3D", "G3DSA"], ["PROSITEPATTERNS", "PS0"], ["SUPERFAMILY", "SSF"],
                  ["CDD", "CD"], ["TIGRFAM", "TIGR"], ["PIRSF", "PIRSF"],
                  ["PRINTS", "PR"], ["HAMAP", "MF"], ["PRODOM", "PD"],
                  ["SFLD", "SFLD"], ["PANTHER", "PTHR"], ["NA", "NAN"]]

    d = {lst[0]: lst[1] for lst in lst_of_lst}
    for k, v in d.items():
        if string.startswith(v):
            return k


def get_cc_len(fn):
    """
    Retrieves the number of nodes contained in a CC
    """

    l, cl = {}, {}

    for k, v in fn.items():
        l[k] = len(v)

    for k, v in l.items():
        if v not in cl.keys():
            cl[v] = 0
        cl[v] += 1

    return cl


def create_list_from_list_of_string(list_of_string):
    """
    Transform a list of large string, into a list of strings
    """
    lst = []
    for i in list_of_string:
        if ',' in i:
            lst.extend(i.split(","))
        else:
            lst.append(i)
    return lst


def standardize(data):
    """
    standardize the data
    """
    try:
        return create_list_from_list_of_string(get_data_from_string(data))
    except:
        return data


def get_data_from_cc(g_cc, columns, out_files, ov, ide, ev):
    """
    Retrieves all wanted data from a connected component (CC)
    """
    d = {}
    for k, cc in enumerate(g_cc):
        cc_n = k + 1
        d[cc_n] = {}
        for j in columns:
            cc.vs[j] = standardize(cc.vs[j])
            d[cc_n][j] = get_percent(cc.vs[j])
    save_data_dict(d, columns, out_files, ov, ide, ev)


def find_output(c, outputs, ov, ide, ev):
    """
    Retrieves the output file name based on column name, overlap, identity, and evalue parameters
    """
    if ov == -1 or ide == -1 or ev == -1:
        for out in outputs:
            if c in out:
                return out

    for out in outputs:
        if c in out and str(ov) in out and str(ide) in out and str(ev) in out:
            return out


def save_data_dict(d, columns, outputs, ov, ide, ev):
    """
    Saves results in output files
    """
    for c in columns:
        o = find_output(c, outputs, ov, ide, ev)
        tmp = []
        for k, v in d.items():
            tmp.extend([str(k), str(k2), str(v2)] for k2, v2 in v[c].items())
        with open(o, "w") as f:
            print(f"CC\t{c}\tPercentage", file=f)
            for i in tmp:
                print("\t".join(i), file=f)


def get_data_from_string(lst_of_lst):
    """
    Retrieves a list of all database ID contained
    in a list of long string of database ID
    """

    a = []

    for lst in lst_of_lst:

        if lst == "nan":
            a.append(l)

        else:
            lst = str(lst)
            sl = lst.split("|")

            a.extend(iter(sl))
    return a


def get_homogeneity_score(d):
    """
    Retrieves a homogeneity score from a dictionary
    """

    le = sum(v for k, v in d.items())

    u = len(d)

    if u == 1:
        return 1

    return round(1 - (u / le), 3)


def get_percent(lst):
    """
    Retrieves a percentage dictionary of all elements in a list
    """
    tot = 0
    d = {}

    for i in lst:
        if type(i) == str:
            k = i.split(",")
            sl = list(k)
            for j in sl:
                if j not in d.keys():
                    d[j] = 0
                d[j] += 1
                tot += 1

        else:
            if i not in d.keys():
                d[i] = 0
            d[i] += 1
            tot += 1

    return {k: percentage(v, tot) for k, v in d.items()}


def get_hom_score_data_from_cc(g_cc, columns, outputs, ov, ide, ev):
    """
    Gets a dictionary filled with homogeneity scores from CC
    """
    d = {}
    for k, cc in enumerate(g_cc):
        cc_n = k + 1
        d[cc_n] = {
            f"homogeneity_score_{j}": get_homogeneity_score(Counter(cc.vs[j]))
            for j in columns
        }

    save_hom_score(d, columns, outputs, ov, ide, ev)


def get_homogeneity_key(v, c):
    """
    Retrieves the key of a dict
    """
    for k in v.keys():
        if c in k:
            return k


def save_hom_score(d, columns, outputs, ov, ide, ev):
    """
    Saves homogeneity scores in output files
    """
    for c in columns:
        tmp_hom = []
        out = find_output(c, outputs, ov, ide, ev)
        for k, v in d.items():
            v_homog = get_homogeneity_key(v, c)
            tmp_hom.append([str(k), str(v[v_homog])])
        with open(out, "w") as f:
            print("CC\tHomogeneity_score", file=f)
            for i in tmp_hom:
                print("\t".join(i), file=f)


def get_set_from_col(c, g_cc):
    """
    Gets a set of column names
    """
    set_l = set()
    for cc in g_cc:
        cc.vs[c] = standardize(cc.vs[c])
        for d_col in cc.vs[c]:
            if d_col not in set_l:
                if c in ["name", "peptides"]:
                    set_l.add(col.split(".")[0])
                else:
                    set_l.add(col)
    return l


def fill_abund_dict(d, g_cc, c):
    """
    Fills the abundance matrix
    """
    for index, cc in enumerate(g_cc):
        cc.vs[c] = standardize(cc.vs[c])
        d[index] = {'CC': index}
        for d_col in cc.vs[c]:
            if c in ["name", "peptides"]:
                if d_col not in d:
                    d[index][d_col.split(".")[0]] = 0
                d[index][d_col.split(".")[0]] += 1
            else:
                if d_col not in d:
                    d[index][d_col] = 0
                d[index][d_col] += 1
    return d


def get_fieldnames(c, g_cc):
    """
    Retrieves the fieldnames of the output file
    """
    tmp = set()
    for cc in g_cc:
        cc.vs[c] = standardize(cc.vs[c])
        for col in cc.vs[c]:
            if col not in tmp:
                if c in ["name", "peptides"]:
                    tmp.add(col.split(".")[0])
                else:
                    tmp.add(col)

    return ['CC', *[i for i in tmp]]


def save_matrix(d, output, fieldnames):
    """
    Saves the abundance matrix in output file
    """
    writer = csv.DictWriter(open(output, "w"), fieldnames=fieldnames)
    writer.writeheader()
    for k, v in d.items():
        writer.writerow(v)


def get_abund_mat(g_cc, columns, outputs, ov, ide, ev):
    """
    Retrieves the abundance matrix and saves it in the output file
    """
    d = {}
    for c in columns:
        output = find_output(c, outputs, ov, ide, ev)
        fieldnames = get_fieldnames(c, g_cc)
        d = fill_abund_dict(d, g_cc, c)
        save_matrix(d, output, fieldnames)


def find_io(overlap, identity, e_val, edges, vertices, outputs):
    """
    Retrieves the input and output files based on the overlap, the identity, and the e value in
    parameters
    """
    in_files, out_files = [], []
    length = len(edges)

    for i in range(length):
        e = edges[i].split("_")
        v = vertices[i].split("_")
        if f"{overlap}" in e[2] and f"{identity}" in e[3] and f"{e_val}" in e[4]:
            in_files.append(edges[i])
        if f"{overlap}" in v[2] and f"{identity}" in v[3] and f"{e_val}" in v[4]:
            in_files.append(vertices[i])

    length = len(outputs)

    for i in range(length):
        f = outputs[i].split("_")
        if f"{overlap}" in f[0] and f"{identity}" in f[1] and f"{e_val}" in f[2]:
            out_files.append(outputs[i])

    return in_files, out_files


def main():
    """
    Main program function
    """

    if not snakemake.params.similarity:
        for ov in snakemake.params.overlap:
            for ide in snakemake.params.identity:
                for ev in snakemake.params.eval:
                    in_files, out_files = find_io(ov, ide, ev, snakemake.input.edges,
                                                  snakemake.input.vertices, snakemake.output.rslts)

                    # creates pandas dataframe of edges and nodes
                    edges = pd.read_csv(in_files[0], sep=";")
                    nodes = pd.read_csv(in_files[1], sep=";", low_memory=False)

                    # create an igraph Graph
                    g = ig.Graph.DataFrame(edges, directed=False)

                    # remove isolated nodes
                    if snakemake.params.isolated == "yes":
                        to_del = [v.index for v in g.vs if v.degree() == 0]
                        g.delete_vertices(to_del)

                    for index, i in enumerate(snakemake.params.columns):
                        t = list(nodes["name"]) if index == 0 else list(nodes[i])
                        g.vs[i] = t

                    # decompose graph into subgraph
                    g_cc = g.decompose(minelements=snakemake.params.neighbours)

                    # get a dictionary containing all data from CC and save it
                    get_data_from_cc(g_cc, snakemake.params.columns, out_files, ov, ide, ev)

                    # get homogeneity score from each CC
                    get_hom_score_data_from_cc(g_cc, snakemake.params.columns,
                                               snakemake.output.homscore, ov, ide, ev)

                    # get an abundance matrix from CC and save it
                    get_abund_mat(g_cc, snakemake.params.columns, snakemake.output.abund_mat,
                                  ov, ide, ev)

                    ig.Graph.write_graphml(g, f"results/ssn_graph_{ov}_{ide}_{ev}")

    else:
        ov, ide, ev = -1, -1, -1
        for f in snakemake.params.similarity:
            g = ig.Graph.Read_GraphML(str(f))
            out_files = [f"results/{f.split('/')[-1]}_{c}_results"
                         for c in snakemake.params.columns]

            hom_score_files = [f"results/{f.split('/')[-1]}_{c}_homogeneity_score"
                               for c in snakemake.params.columns]

            abund_mat_files = [f"results/{f.split('/')[-1]}_{c}_abundance_matrix"
                               for c in snakemake.params.columns]

            if snakemake.params.isolated == "yes":
                to_del = [v.index for v in g.vs if v.degree() == 0]
                g.delete_vertices(to_del)

            # decompose graph into subgraph
            g_cc = g.decompose(minelements=snakemake.params.neighbours)

            # get a dictionary containing all data from CC
            get_data_from_cc(g_cc, snakemake.params.columns, out_files, ov, ide, ev)

            # get homogeneity score from each CC
            get_hom_score_data_from_cc(g_cc, snakemake.params.columns,
                                       hom_score_files, ov, ide, ev)

            # get an abundance matrix from CC and save it
            get_abund_mat(g_cc, snakemake.params.columns, abund_mat_files,
                          ov, ide, ev)


if __name__ == '__main__':
    main()

# def get_entropy(d):
#    """
#    Retrieves the Shannon entropy of dictionary values
#
#    Parameters
#    ----------
#
#        d : a dictionary containing values on databases
#
#    Returns
#    -------
#
#        en : the entropy of the dict values
#    """
#    n = len(d)
#    en = 0
#    for v in d.values():
#        val = v / 100
#        if val == 1:
#            return val
#        elif val == 0:
#            return 0
#        else:
#            en += -val * math.log(val, n)
#    return round(en, 3)
