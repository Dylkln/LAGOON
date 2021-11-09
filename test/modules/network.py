import igraph as ig
import argparse
import pandas as pd
import math
import os


def arguments():
    """
    set arguments
    """

    parser = argparse.ArgumentParser()

    # Mandatory arguments
    parser.add_argument("-e", "--edges_file", dest="edges_file",
                        type=str, required=True, nargs="+",
                        help="")
    parser.add_argument("-v", "--vertices_file", dest="vertices_file",
                        type=str, required=True, nargs="+",
                        help="")
    parser.add_argument("-n", "--neighbours", dest="neighbours",
                        type=int, required=False, default=3,
                        help="")
    parser.add_argument("-sn", "--similarity_network", dest="similarity_network",
                        type=str, required=False,
                        help="")
    parser.add_argument("-i", "--isolated", dest="isolated",
                        type=str, required=False, default="No",
                        help="")
    parser.add_argument("-cn", "--columns", dest="columns",
                        type=str, required=False, default=None, nargs="+",
                        help="the columns names used in attributes")

    return parser.parse_args()


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


def get_data_from_cc(g_cc, columns):
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

    return d

def save_data_dict(d, columns, cov, ident):
    for c in columns:
        f = open(f"results/pcov{cov}_pident{ident}_ssn_{c}_results", "w")
        for k, v in d.items():
            for k2, v2 in v[c].items():
                f.write(f"{k}\t{k2}\t{v2}\n")



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


def main():
    """
    Main program function
    """

    # get arguments
    args = arguments()

    if not args.similarity_network:

        # get number of files
        l_e = len(args.edges_file)

        for i in range(l_e):

            # get coverage percentage
            cov = args.edges_file[i].split("_")[2].strip("pcov")

            # get identity percentage
            ident = args.edges_file[i].split("_")[3].split(".")[0].strip("pident")

            # creates pandas dataframe of edges and nodes
            edges = pd.read_csv(args.edges_file[i], sep=";")
            nodes = pd.read_csv(args.vertices_file[i], sep=";", low_memory=False)

            # create an igraph Graph
            g = ig.Graph.DataFrame(edges, directed=False)

            # remove isolated nodes
            if args.isolated == "yes":
                to_del = [v.index for v in g.vs if v.degree() == 0]
                g.delete_vertices(to_del)


            for index, i in enumerate(args.columns):
                t = [j for j in nodes["name"]] if index == 0 else [j for j in nodes[i]]
                g.vs[i] = t
            # decompose graph into subgraph
            g_cc = g.decompose(minelements=args.neighbours)

            # get number of connected components
    #        nb_of_subgraph = len(g_cc)

            # get a dictionary containing all data from CC
            save_data_dict(get_data_from_cc(g_cc, args.columns), args.columns, cov, ident)

            # get trophy count, unique trophy, homogeneity and entropy
            tr_c, u_tr, hi, f_en = get_data_from_cc_dicts(tr, db_sp, fn_p)

            # get the number of IDs contained in each connected component
            cl = get_cc_len(fn)

            # get abundance distribution
            ab_d = get_abund_distrib(ab)

            # check path
            p = f"../results/{cov}_{ident}"

            if not os.path.exists(p):
                os.mkdir(p)

            # save all values extracted from the graph
            output = f"../results/{cov}_{ident}/abund_matrix_{cov}_{ident}.tsv"
            save_dict_of_dict(ab, output)

            output = f"../results/{cov}_{ident}/abund_matrix_distrib_{cov}_{ident}.tsv"
            save_dict(ab_d, output)

            output = f"../results/{cov}_{ident}/function_percentage_{cov}_{ident}.tsv"
            save_dict_of_dict(fn_p, output)

            output = f"../results/{cov}_{ident}/phylum_percentage_{cov}_{ident}.tsv"
            save_dict_of_dict(tp_p, output)

            output = f"../results/{cov}_{ident}/genus_percentage_{cov}_{ident}.tsv"
            save_dict_of_dict(tg_p, output)

            output = f"../results/{cov}_{ident}/trophy_percentage_{cov}_{ident}.tsv"
            save_dict_of_dict(tr_p, output)

            output = f"../results/{cov}_{ident}/cc_nb_{cov}_{ident}.txt"
            save_dict(cl, output)

            output = f"../results/{cov}_{ident}/trophy_count_{cov}_{ident}.tsv"
            save_dict_of_dict(tr_c, output)

            output = f"../results/{cov}_{ident}/unique_trophy_{cov}_{ident}.tsv"
            save_dict(u_tr, output)

            output = f"../results/{cov}_{ident}/homogeneity_index_{cov}_{ident}.tsv"
            save_dict_of_dict(hi, output)

            #        output = f"../results/{cov}_{ident}/entropy_by_function_{cov}_{ident}.tsv"
            #        save_dict_with_cc(f_en, output)

            # save the graph into a graphml format
            g.write(f=f"../results/{cov}_{ident}/graph_ssn_{cov}_{ident}", format="graphml")


if __name__ == '__main__':
    main()
