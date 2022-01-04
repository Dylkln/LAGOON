import math

import igraph as ig
import pandas as pd


def percentage(part, whole):
    percent = 100 * float(part) / float(whole)
    return round(percent, 2)


def db(string):
    """
    Retrieves the name of the database based on the string provided in argument

    Parameters
    ----------

        string : a database ID

    Returns
    -------

        k : the database associated to the database ID
    """
    ll = [["PFAM", "PF"], ["SMART", "SM"], ["PROSITEPROFILES", "PS5"],
          ["GENE3D", "G3DSA"], ["PROSITEPATTERNS", "PS0"], ["SUPERFAMILY", "SSF"],
          ["CDD", "cd"], ["TIGRFAM", "TIGR"], ["PIRSF", "PIRSF"],
          ["PRINTS", "PR"], ["HAMAP", "MF"], ["PRODOM", "PD"],
          ["SFLD", "SFLD"], ["PANTHER", "PTHR"], ["NA", "nan"]]

    d = {l[0]: l[1] for l in ll}
    for k, v in d.items():
        if string.startswith(v):
            return k


def get_abund(l):
    """
    Retrieves the abundance of each element of a list

    Parameters
    ----------

        l : a list of transcript

    Returns
    -------

        d : a dictionary representing the abundance of each trancript in the list
    """
    d = {}

    for i in l:
        iid = "-".join(i.split("-")[0:2])

        if iid not in d.keys():
            d[iid] = 0
        d[iid] += 1

    return d


def get_abund_distrib(d):
    """
    Retrieves the distribution of the abundance of a dictionary

    Parameters
    ----------

        d : a dictionary with the abundance of each transcript

    Returns
    -------

        ad : a dictionary containing the distribution of the abundance
    """

    ad = {}

    for k, v in d.items():
        if len(v) == 1:

            for k2 in v.keys():
                if k2 not in ad.keys():
                    ad[k2] = 0
                ad[k2] += 1
    return ad


def get_cc_len(fn):
    """
    Retrieves the number of nodes contained in a CC

    Parameters
    ----------

        fn : a dictionary containing all Databases in the CC

    Returns
    -------

        cl : a dictionary containing all nodes number of each CC
    """

    l, cl = {}, {}

    for k, v in fn.items():
        l[k] = len(v)

    for k, v in l.items():
        if v not in cl.keys():
            cl[v] = 0
        cl[v] += 1

    return cl


def get_count(l):
    """
    Retrieves the count of each value on a list

    Parameters
    ----------

        l : a list

    Returns
    -------

        c : a dictionary containing the count of each value of the list
    """

    c = {}

    for i in l:
        if i not in c.keys():
            c[i] = 0
        c[i] += 1

    return c


def get_data_from_cc(g_cc, columns, out_files):
    """
    Retrives all wanted data from a connected component (CC)

    Parameters
    ----------

        g_cc: a connected componant of a graph

    Returns
    -------

        fn : a dictionary containing all Databases in the CC
        tp : a dictionary containing all Phylums in the CC
        tg : a dictionary containing all Genus in the CC
        tr : a dictionary containing all Trophies information in the CC
        ab : a dictionary containing the abundance of transcript in the CC
    """
    d = {}
    for k, cc in enumerate(g_cc):
        cc_n = k + 1
        d[cc_n] = {}
        for i, j in enumerate(columns):
            if "database" in j:
                f = get_db_id(cc.vs[j])
                d[cc_n][j] = get_percent(get_db(f))
            d[cc_n][j] = get_percent(cc.vs[j])
    save_data_dict(d, columns, out_files)


def find_output(c, outputs):
    for out in outputs:
        if c in out:
            return out


def save_data_dict(d, columns, outputs):
    for c in columns:
        o = find_output(c, outputs)
        tmp = []
        for k, v in d.items():
            for k2, v2 in v[c].items():
                tmp.append([str(k), str(k2), str(v2)])
        with open(o, "w") as f:
            print(f"CC\t{c}\tPercentage", file=f)
            for i in tmp:
                print("\t".join(i), file=f)



def get_data_from_cc_dicts(tr, db_sp, fn_p):
    """
    Retrieves different data from the several connected components dictionaries

    Parameters
    ----------

        tr : a dictionary containing all Trophies information in the CC
        db_sp :
        fn_p : a dictionary containing the percentage of all Databases
        in the CC

    Returns
    -------

        tr_c : a dictionary containing the number of trophies in each CC
        u_tr : a dictionary containing the number CC with only one trophy
        hi : a dictionary containing the homogeneity score from each Database
        of each CC
        f_en : a dictionary containing the entropy of all Databases on each CC

    """

    tr_c, u_tr, hi, f_en = {}, {}, {}, {}

    for k, v in tr.items():
        if k not in tr_c.keys():
            tr_c[k] = get_count(v)

    for k, v in tr_c.items():
        c = len(v)
        if c == 1:
            for key in v.keys():
                if key not in u_tr.keys():
                    u_tr[key] = 0
                u_tr[key] += 1
        else:
            if "not_unique" not in u_tr.keys():
                u_tr["not_unique"] = 0
            u_tr["not_unique"] += 1

    for k, v in db_sp.items():
        for k2, v2 in v.items():
            if k not in hi.keys():
                hi[k] = {}
            if k2 not in hi[k].keys():
                hi[k][k2] = get_homogeneity_score(v2)

        for ke, va in fn_p.items():
            if ke not in f_en.keys():
                f_en[ke] = get_entropy(va)

    return tr_c, u_tr, hi, f_en


def get_data_percent(fn, tp, tg, tr):
    """
    Retrieves data percentage of each information extracted from a CC

    Parameters
    ----------

        fn : a dictionary containing all Databases in the CC
        tp : a dictionary containing all Phylums in the CC
        tg : a dictionary containing all Genus in the CC
        tr : a dictionary containing all Trophies information in the CC

    Returns
    -------

        fn_p : a dictionary containing the percentage of all Databases
        in the CC
        tp_p : a dictionary containing the percentage of all Phylums
        in the CC
        tg_p : a dictionary containing the percentage of all Genus
        in the CC
        tr_p : a dictionary containing the percentage of all Trophies
        information in the CC
    """

    fn_p, tp_p, tg_p, tr_p = {}, {}, {}, {}

    for k in fn.keys():
        fn_p[k] = get_percent(fn[k])
        tp_p[k] = get_percent(tp[k])
        tg_p[k] = get_percent(tg[k])
        tr_p[k] = get_percent(tr[k])

    return fn_p, tp_p, tg_p, tr_p


def get_db(l):
    """
    Retrives a list of all database contained in a CC

    Parameters
    ----------

        l : a list of database IDs

    Returns
    -------

        db_l : a list of all database contained in a CC
    """

    return [db(i) for i in l]


def get_db_id(l_o_l):
    """
    Retrieves a list of all database ID contained
    in a list of long string of database ID

    Parameters
    ----------

        l_o_l : a list of long string of database ID

    Returns
    -------

        a : a list of all database ID contained in l_o_l
    """

    a = []

    for l in l_o_l:

        if l == "nan":
            a.append(l)

        else:
            l = str(l)
            sl = l.split("|")

            for i in sl:
                a.append(i)

    return a


def get_db_sp(l):
    """
    Retrieves a dictionary containing all annotations for each CC

    Parameters
    ----------

        l : a list of IDs

    Returns
    -------

        db_sp : a dictionary containing all IDs by database
    """

    db_sp = {}
    lt = transform_list(l)

    for i in lt:
        name = db(i)
        if name not in db_sp.keys():
            db_sp[name] = {}
        if i not in db_sp[name].keys():
            db_sp[name][i] = 0
        db_sp[name][i] += 1

    return db_sp


def get_entropy(d):
    """
    Retrieves the Shannon entropy of dictionary values

    Parameters
    ----------

        d : a dictionary containing values on databases

    Returns
    -------

        en : the entropy of the dict values
    """
    n = len(d)
    en = 0
    for v in d.values():
        val = v / 100
        if val == 1:
            return val
        elif val == 0:
            return 0
        else:
            en += -val * math.log(val, n)
    return round(en, 3)


def get_homogeneity_score(d):
    """
    Retrieves an homogeneity score from a dictionary

    Parameters
    ----------

        d : a dictionary containing values on databases

    Returns
    -------

        hi : a dictionary containing the homogeneity score for each database
    """

    le = sum(v for k, v in d.items())

    u = len(d)

    if u == 1:
        return 1

    return round(1 - (u / le), 3)


def get_percent(l):
    """
    Retrieves a percentage dictionary of all elements in a list

    Parameters
    ----------

        l : a list of elements

    Returns
    -------

        d_percent : a dictionary containing all percentages of each element of
        the list given in argument
    """
    tot = 0
    d = {}

    for i in l:
        if type(i) == str:
            k = i.split(",")
            sl = [j for j in k]
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


def save_dict_with_cc(d, output):
    """
    Save all data contained in a dictionary in a file

    Parameters
    ----------

        d : a dictionary
        output : the file where to save the data from the dictionary
    """

    with open(output, "w") as f:
        for k, v in d.items():
            f.write(f"CC_{k}\t{v}\n")


def save_dict(d, output):
    """
    Save all data contained in a dictionary in a file

    Parameters
    ----------

        d : a dictionary
        output : the file where to save the data from the dictionary
    """

    with open(output, "w") as f:
        for k, v in d.items():
            f.write(f"{k}\t{v}\n")


def save_dict_of_dict(d, output):
    """
    Save all data contained in a dictionary of dictionary in a file

    Parameters
    ----------

        d : a dictionary of dictionary
        output : the file where to save the data from the dictionary
    """

    with open(output, "w") as f:
        for k, v in d.items():
            for k2, v2 in v.items():
                f.write(f"CC_{k}\t{k2}\t{v2}\n")


def save_list_in_dict(d, output):
    """
    Save all data contained in a dictionary of lists in a file

    Parameters
    ----------

        d : a dictionary of lists
        output : the file where to save the data from the dictionary
    """

    with open(output, "w") as f:
        for k, v in d.items():
            for i in v:
                f.write(f"CC_{k}\t{i}\n")


def save_value(v, output):
    """
    Save data from a variable in a file

    Parameters
    ----------

        v : a variable containg a value
        output : the file where to save the data from the variable
    """
    with open(output, "w") as f:
        f.write(v)


def transform_list(l):
    """
    Retrieves a list containing all values of another list separated by "|"

    Parameters
    ----------

        l : a list of long string separated by "|"

    Returns
    -------

        lt : list of all values contained in l
    """

    lt = []
    for i in l:
        j = i.split("|")
        for k in j:
            lt.append(k)
    return lt


def find_io(overlap, identity, eval, edges, vertices, outputs):
    in_files = []
    length = len(edges)

    for i in range(length):
        if f"{overlap}" in edges[i] and f"{identity}" in edges[i] and f"{eval}" in edges[i]:
            in_files.append(edges[i])
        if f"{overlap}" in vertices[i] and f"{identity}" in vertices[i] and f"{eval}" in \
                vertices[i]:
            in_files.append(vertices[i])

    out_files = [
        f
        for f in outputs
        if f"{overlap}" in f and f"{identity}" in f and f"{eval}" in f
    ]

    return in_files, out_files


def main():
    """
    Main program function
    """

    similarity = snakemake.params.similarity
    if not similarity:

        in_edges = snakemake.input.edges
        in_vertices = snakemake.input.vertices
        outputs = snakemake.output
        overlap = snakemake.params.overlap
        identity = snakemake.params.identity
        neighbours = snakemake.params.neighbours
        columns = snakemake.params.columns
        isolated = snakemake.params.isolated
        eval = snakemake.params.eval

        for ov in overlap:
            for id in identity:
                for ev in eval:
                    in_files, out_files = find_io(ov, id, ev, in_edges, in_vertices, outputs)

                    # creates pandas dataframe of edges and nodes
                    edges = pd.read_csv(in_files[0], sep=";")
                    nodes = pd.read_csv(in_files[1], sep=";", low_memory=False)

                    # create an igraph Graph
                    g = ig.Graph.DataFrame(edges, directed=False)

                    # remove isolated nodes
                    if isolated == "yes":
                        to_del = [v.index for v in g.vs if v.degree() == 0]
                        g.delete_vertices(to_del)


                    for index, i in enumerate(columns):
                        t = [j for j in nodes["name"]] if index == 0 else [j for j in nodes[i]]
                        g.vs[i] = t
                    # decompose graph into subgraph
                    g_cc = g.decompose(minelements=neighbours)

                    # get number of connected components
            #        nb_of_subgraph = len(g_cc)

                    # get a dictionary containing all data from CC
                    get_data_from_cc(g_cc, columns, out_files)

if __name__ == '__main__':
    main()
