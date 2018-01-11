"""
Routines to read and write lattices from and to files.
"""

import numpy as np

from lattice import Lattice, Site

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
        wndf.write("0")
