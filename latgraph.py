"""
Main program. Convert lattice file formats / relabel adjacency graphs / generate lattices.
"""

import sys
import argparse
import textwrap

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import fileio
import tubegen


# List all available generators here.
# Each generator must have a 'run' function that takes command line arguments
# as its sole parameter and returns a lattice.
GENERATORS = {"tube": tubegen}

def define_parser():
    "Define the main command line argument parser"

    parser = argparse.ArgumentParser(prog="latgraph",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""\
            Manipulate, generate convert lattices.

            Can either load a lattice from file (-i) or generate a new lattice (-g).
                - Use -i filename to load load a file with a lattice.
                  The type is deduced from the extension

                - Use -g generator to run a lattice generator. All arguments you
                  specify after that are passed to the generator.
                  You can get generator specific help via -g generator -h
                  Available generators are:
                      {}
                                     """.format(", ".join(GENERATORS.keys()))),
                                     epilog="See https://github.com/jl-wynen/latgraph")

    cmd_group = parser.add_mutually_exclusive_group(required=True)
    cmd_group.add_argument("-i", "--input", metavar="infile", help="Input lattice file")
    cmd_group.add_argument("-g", "--generate", metavar="generator",
                           nargs=argparse.REMAINDER,
                           help="Generate a lattice using the specified generator")

    parser.add_argument("-o", "--output", metavar="outfile", help="Output lattice file")
    parser.add_argument("-l", "--labels", metavar="labelfile",
                        help="Relabel the graph according to embedding in given file")
    parser.add_argument("-m", "--method", help="Method for relabelling",
                        default="anticlockwise", metavar="method")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="Plot the graph (after relabelling)")
    parser.add_argument("-P", "--plot-labels", action="store_true",
                        help="Plot the label graph")
    parser.add_argument("-a", "--plot_adjacency", action="store_true",
                        help="Show the adjacency graph")

    return parser

def parse_args():
    "Parse command line arguments."
    parser = define_parser()
    args = parser.parse_args()

    if args.generate:
        gen = args.generate[0].lower()
        if not gen in GENERATORS:
            parser.error("Unknown generator: {}".format(gen))
        args.generator = gen
        args.generator_args = args.generate[1:]

    return args

def show_lattice(lat, labels=None):
    "Show the lattice in its own figure using its 3D or 2D embedding."

    fig = plt.figure(figsize=(10, 10))
    if len(lat.sites[0].pos) == 3:
        ax = fig.add_subplot(111, projection="3d")
    else:
        ax = fig.add_subplot(111)
    ax.set_title(lat.name)
    ax.axis("equal")

    for site in lat:
        for neigh in site.neighbours:
            if neigh > site.idx:
                ax.plot(*zip(site.pos, lat.sites[neigh].pos), c="C0")
        if labels is not None:
            ax.text(*site.pos, labels[site.idx])

    centre = lat.centre()
    ax.scatter((centre[0], ), (centre[1], ), marker="x", c="k")
    fig.tight_layout()

def show_adjacency_matrix(lat, name):
    "Show the adjacency matrix in its own figure."

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    ax.set_title(name)
    ax.axis("equal")

    # fill matrix with hopping strengths
    adj = np.zeros((len(lat), len(lat)), dtype=int)
    for i, site in enumerate(lat):
        for j, hop in zip(site.neighbours, site.hopping):
            adj[i, j] = hop

    ax.matshow(adj)
    fig.tight_layout()

def load_lattice(fname):
    lat = fileio.read(fname)
    if not lat.check_consistency():
        print("Input lattice is inconsistent")
        sys.exit(1)
    if not lat.name:
        lat.name = fname
    return lat

def generate_lattice(gen, args):
    return GENERATORS[gen].run(args)

def main():
    args = parse_args()

    if args.input:
        lat = load_lattice(args.input)
    else:
        lat = generate_lattice(args.generator, args.generator_args)

    if args.labels:
        labelLat = fileio.read(args.labels)
        if not labelLat.check_consistency():
            print("Label lattice is inconsistent")
            sys.exit(1)
        labels = labelLat.label_graph(args.method)

        if args.plot_labels:
            show_lattice(labelLat, "Label graph ({})".format(args.labels), labels)

        lat.relabel(labels)

    if args.plot:
        if args.labels:
            show_lattice(lat, labels)
        else:
            show_lattice(lat, np.arange(len(lat)))

    if args.plot_adjacency:
        show_adjacency_matrix(lat, "Adjacency matrix ({})".format(args.infile))

    if args.output:
        fileio.write(args.output, lat)

    if args.plot or args.plot_labels or args.plot_adjacency:
        plt.show()

if __name__ == "__main__":
    main()
