"""
Routines to read and write lattices from and to files.
"""

from pathlib import Path

import numpy as np
import yaml

from lattice import Lattice, Site

# ---------------------- writegraphnd ----------------------

def _wnd_to_lattice(graph):
    "Turn a graph read from a w3d or w2d file into a Lattice object."
    lat = Lattice()
    lat.sites = [Site(idx, pos, neigh, [1]*len(neigh)) for idx, pos, neigh in graph]
    return lat

def read_w3d(fname):
    "Read a Lattice from a writegraph3d file."

    graph = []
    with open(fname, "r") as w3df:
        # check header
        header = w3df.readline().strip()
        if header != ">>writegraph3d<<":
            raise RuntimeError(("Header of file '{}' not correct. Got '{}',"
                                + " expected '>>writegraph3d<<'").format(fname, header))

        # read data
        for line in w3df.readlines():
            if line.strip() == "0":
                break

            numbers = [s for s in line.split() if s]
            graph.append((int(numbers[0])-1,
                          np.array((float(numbers[1]), float(numbers[2]), float(numbers[3]))),
                          list(map(lambda x: int(x)-1, numbers[4:]))))

    return _wnd_to_lattice(graph)

def read_w2d(fname):
    "Read a Lattice from a writegraph2d file."

    graph = []
    with open(fname, "r") as w2df:
        # check header
        header = w2df.readline().strip()
        if header != ">>writegraph2d<<":
            raise RuntimeError(("Header of file '{}' not correct. Got '{}',"
                                + " expected '>>writegraph2d<<'").format(fname, header))

        # read data
        for line in w2df.readlines():
            if line.strip() == "0":
                break

            numbers = [s for s in line.split() if s]
            graph.append((int(numbers[0])-1,
                          np.array((float(numbers[1]), float(numbers[2]))),
                          list(map(lambda x: int(x)-1, numbers[3:]))))

    return _wnd_to_lattice(graph)

def write_wnd(fname, lat):
    "Write a 2D / 3D lattice to a writegraph2d / writegraph3d file."

    with open(fname, "w") as wndf:
        wndf.write(">>writegraph{}d<<\n".format(len(lat.sites[0].pos)))

        for site in lat:
            wndf.write("{:4d}\t\t{:s}\t\t{:s}\n" \
                       .format(site.idx+1,
                               "\t".join(map(str, site.pos)),
                               " ".join(map(lambda x: str(x+1), site.neighbours))))
        wndf.write("0\n")


# ---------------------- yaml ----------------------


def read_yaml(fname):
    "Read a Lattice from a YAML file."
    with open(fname, "r") as yamlf:
        return yaml.safe_load(yamlf)

def _append_adj_hop(site, neigh, hop):
    "Append a connection to the neighbour and hopping lists."
    if neigh not in site.neighbours:
        site.neighbours.append(neigh)
        site.hopping.append(hop)

def _parse_yaml(nt, adjacency, hopping, positions, name="", comment=""):
    "Parse a !lattice YAML node."

    # turn hopping into a list if it isn't already
    if not isinstance(hopping, list):
        hopping = [hopping]*len(adjacency)
    elif len(adjacency) != len(hopping):
        raise RuntimeError("Lengths of adjacency matrix and list of hopping strengths do not match")

    # create the lattice with partly empty sites
    lat = Lattice(name=name, comment=comment, nt=nt)
    lat.sites = [Site(i, np.array(pos), [], []) for i, pos in enumerate(positions)]

    # fill in neighour and hopping lists
    for (i, j), hop in zip(adjacency, hopping):
        if i > len(hopping) or j > len(hopping):
            raise RuntimeError("Site index out of range: [{}, {}]".format(i, j))

        _append_adj_hop(lat.sites[i], j, hop)
        _append_adj_hop(lat.sites[j], i, hop)

    return lat

yaml.add_constructor("!lattice",
                     lambda loader, node: \
                     _parse_yaml(**loader.construct_mapping(node, deep=True)),
                     Loader=yaml.SafeLoader)


def write_yaml(fname, lat):
    "Write a Lattice to a YAML file."
    with open(fname, "w") as yamlf:
        return yaml.dump(lat, yamlf)

def _yaml_represent_lattice(dumper, lat):
    "Create a YAML representation of a Lattice using a !lattice node."
    adj, hopping = zip(*[([site.idx, neigh], hop) for site in lat
                         for (neigh, hop) in zip(site.neighbours, site.hopping)
                         if neigh > site.idx])
    pos = [list(map(float, site.pos)) for site in lat]
    return dumper.represent_mapping("!lattice",
                                    {"name": lat.name,
                                     "comment": lat.comment,
                                     "nt": lat.nt,
                                     "adjacency": list(adj),
                                     "hopping": list(hopping),
                                     "positions": pos},
                                    flow_style=False)

yaml.add_representer(Lattice, _yaml_represent_lattice)


# ---------------------- wrappers ----------------------


def read(fname):
    "Read a lattice from a file. The format is deduced from the file extension."

    suffix = Path(fname).suffix
    if suffix == ".yaml" or suffix == ".yml":
        return read_yaml(fname)
    if suffix == ".w3d":
        return read_w3d(fname)
    if suffix == ".w2d":
        return read_w2d(fname)
    raise ValueError("Unknown file extension: '{}'".format(suffix))

def write(fname, lat):
    "Write a lattice to file. The format is deduced from the file extension."

    suffix = Path(fname).suffix
    if suffix == ".yaml" or suffix == ".yml":
        return write_yaml(fname, lat)
    if suffix == ".w3d" or suffix == ".w2d":
        return write_wnd(fname, lat)
    raise ValueError("Unknown file extension: '{}'".format(suffix))
