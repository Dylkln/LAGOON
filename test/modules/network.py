from collections import Counter
import igraph as ig
import pandas as pd


def percentage(part, whole):
    """
    Retrieves the percentage
    """
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
          ["CDD", "CD"], ["TIGRFAM", "TIGR"], ["PIRSF", "PIRSF"],
          ["PRINTS", "PR"], ["HAMAP", "MF"], ["PRODOM", "PD"],
          ["SFLD", "SFLD"], ["PANTHER", "PTHR"], ["NA", "NAN"]]

    d = {l[0]: l[1] for l in ll}
    for k, v in d.items():
        if string.startswith(v):
            return k


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


def create_list_from_list_of_string(list_of_string):
    """
    Transform a list of large string, into a list of strings
    """
    l = []
    for i in list_of_string:
        if ',' in i:
            l.extend(i.split(","))
        else:
            l.append(i)
    return l


def standardize(data):

    try:
        return create_list_from_list_of_string(get_data_from_string(data))
    except:
        return data



def get_data_from_cc(g_cc, columns, out_files):
    """
    Retrieves all wanted data from a connected component (CC)

    Parameters
    ----------
        g_cc : graph connected components
        columns : column names
        out_files : a list of output files
    """
    d = {}
    for k, cc in enumerate(g_cc):
        cc_n = k + 1
        d[cc_n] = {}
        for j in columns:
            cc.vs[j] = standardize(cc.vs[j])
            d[cc_n][j] = get_percent(cc.vs[j])
            d[cc_n][f"homogeneity_score_{j}"] = get_homogeneity_score(Counter(cc.vs[j]))
    save_data_dict(d, columns, out_files)


def find_output(c, outputs):
    for out in outputs:
        if c in out:
            return out


def get_homogeneity_keys(v):
    return [k for k in v.keys() if "homogeneity_score" in k]


def save_data_dict(d, columns, outputs):
    d_h = {}
    os_h = []
    for c in columns:
        o = find_output(c, outputs)
        tmp = []
        for k, v in d.items():
            if "homogeneity_score" in v:
                v_homog = get_homogeneity_keys(v)
                for i in v_homog:
                    d_h[k] = v[i]
            tmp.extend([str(k), str(k2), str(v2)] for k2, v2 in v[c].items())
        with open(o, "w") as f:
            print(f"CC\t{c}\tPercentage", file=f)
            for i in tmp:
                print("\t".join(i), file=f)

    o_hs = "_".join(outputs[0].split("_")[:3]) + "_homogeneity_score"
    with open(o_hs, "w") as f:
        print("CC\tHomogeneity_score", file=f)
        for k, v in d_h.items():
            print(f"{k}\t{v}", file=f)


def get_data_from_string(l_o_l):
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

            a.extend(iter(sl))
    return a


def get_homogeneity_score(d):
    """
    Retrieves a homogeneity score from a dictionary

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


def find_io(overlap, identity, e_val, edges, vertices, outputs):
    """
    Retrieves the input and output files based on the overlap, the identity, and the e value given

    Parameters
    ----------

        overlap : Overlap percentage
        identity : Identity percentage
        e_val : e value
        edges : edges files
        vertices : vertices files
        outputs : a list of all outputs

    returns
    -------

        in_files : input files
        out_files : output files
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
                                                  snakemake.input.vertices, snakemake.output)

                    # creates pandas dataframe of edges and nodes
                    edges = pd.read_csv(in_files[0], sep=";")
                    nodes = pd.read_csv(in_files[1], sep=";", low_memory=False)
                    print(nodes)

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

                    # get a dictionary containing all data from CC
                    get_data_from_cc(g_cc, snakemake.params.columns, out_files)

                    ig.Graph.write_graphml(g, f"results/ssn_graph_{ov}_{ide}_{ev}")

    else:
        for f in snakemake.params.similarity:
            g = ig.Graph.Read_GraphML(str(f))
            o = "_".join(f.split("_")[2:])
            out_files = [
                f'results/{o}' + f"_ssn_{i}_results"
                for i in snakemake.params.columns
            ]


            if snakemake.params.isolated == "yes":
                to_del = [v.index for v in g.vs if v.degree() == 0]
                g.delete_vertices(to_del)

            # decompose graph into subgraph
            g_cc = g.decompose(minelements=snakemake.params.neighbours)

            # get a dictionary containing all data from CC
            get_data_from_cc(g_cc, snakemake.params.columns, out_files)


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
