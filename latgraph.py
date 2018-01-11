"""

"""

import sys
from argparse import ArgumentParser as AParser

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# import lattice
import fileio

def show_lattice(lat, name, labels=None):
    fig = plt.figure(figsize=(10, 10))
    if len(lat.sites[0].pos) == 3:
        ax = fig.add_subplot(111, projection="3d")
    else:
        ax = fig.add_subplot(111)
    ax.set_title(name)

    for site in lat:
        for neigh in site.neighbours:
            if neigh > site.idx:
                ax.plot(*zip(site.pos, lat.sites[neigh].pos), c="C0")
        if labels is not None:
            ax.text(*site.pos, labels[site.idx])

    centre = lat.centre()
    ax.scatter((centre[0], ), (centre[1], ), marker="x", c="k")
    fig.tight_layout()

def parse_args():
    "Parse command line arguments."

    parser = AParser(prog="latgraph",
                     description="""
                     Convert lattice geometries between different file formats.
                     """,
                     epilog="")
    parser.add_argument("input", help="Input lattice file")
    parser.add_argument("-o", "--output", metavar="outfile", help="Output lattice file")
    parser.add_argument("-l", "--labels", metavar="labelfile",
                        help="Relabel the graph according to embedding in given file")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="Plot the graph (after relabelling)")
    parser.add_argument("-P", "--plot-labels", action="store_true",
                        help="Plot the label graph")
    return parser.parse_args()

def main():
    args = parse_args()

    lat = fileio.read(args.input)
    if not lat.check_consistency():
        print("Input lattice is inconsistent")
        sys.exit(1)

    if args.labels:
        labelLat = fileio.read(args.labels)
        if not labelLat.check_consistency():
            print("Label lattice is inconsistent")
            sys.exit(1)
        labels = labelLat.label_graph()

        if args.plot_labels:
            show_lattice(labelLat, "Label graph ({})".format(args.labels), labels)

        lat.relabel(labels)

    if args.plot:
        if args.labels:
            show_lattice(lat, args.input, labels)
        else:
            show_lattice(lat, args.input, np.arange(len(lat)))

    if args.output:
        fileio.write(args.output, lat)

    if args.plot or args.plot_labels:
        plt.show()

if __name__ == "__main__":
    main()
