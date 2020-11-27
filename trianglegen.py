"""
Generator for triangle tilings.
"""

r"""
Generator for triangle tilings.

Generate a trinangle tiling.
The tiling is classified in terms of 'rows' and 'columns' like:
      - columns -
  |    /\/\/\/\
rows   \/\/\/\/
  |    /\/\/\/\
This lattice has 3 rows and 7 columns.

Rows can be stacked in two different ways:
- first = 'up'     /\/\
                   \/\/

- first = 'down'   \/\/\
                   /\/\/

The left edge determines the stycking type.

Boundaries are always open.
"""


import argparse
from enum import Enum, auto
import textwrap

import numpy as np

from lattice import Lattice, Site


def _define_parser():
    "Define command line argument parser for TriangleGen."

    parser = argparse.ArgumentParser(prog="latgraph --generate triangle",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""\
                Generate a triangle tiling.
                                     """))
    parser.add_argument("cols", type=int, help="Number of triangles per row.")
    parser.add_argument("rows", type=int, help="Number of stacked rows.")
    parser.add_argument("first", type=str, choices=("up", "down"),
                        help="Direction the first triangle points in.")
    parser.add_argument("--spacing", type=float, default=1.0, help="Lattice spacing.")
    parser.add_argument("--name", default="", help="Name for the lattice")
    parser.add_argument("--comment", default="", help="Comment on the lattice")
    return parser

def run(in_args):
    "Run TriangleGen from command line arguments."

    args = _define_parser().parse_args(in_args)
    if args.first == "down":
        print("first='down' is not implemented")
        import sys
        sys.exit(1)
    return make_triangle_tiling(args.cols, args.rows, args.first, args.spacing,
                                args.name, args.comment)

def make_triangle_tiling(cols, rows, first, spacing, name, comment):
    assert first == "up"
    lat = Lattice(name, comment)

    def boundary_filter(neigh, row):
        return neigh[0] >= 0 and neigh[0] < _n_vert_cols(neigh[1], cols, first) and \
               neigh[1] >= 0 and neigh[1] <= rows

    for y in range(rows + 1):
        for x in range(_n_vert_cols(y, cols, first)):
            pos = _position(x, y, first, spacing)
            neighbours = tuple(map(lambda t: _total_index(*t, cols, first),
                                   filter(lambda t: boundary_filter(t, y),
                                          _nearest_neighbours(x, y, first))))
            hoppings = [1]*len(neighbours)
            lat.sites.append(Site(_total_index(x, y, cols, first), pos, neighbours, hoppings))

    assert lat.check_consistency()
    return lat

def _ceildiv(a, b):
    return -(-a // b)

def _n_vert_cols(y, cols, first):
    if (first == "up" and y % 2 == 0) or (first == "down" and y % 2 == 1):
        return _ceildiv(cols + 1, 2)
    return (cols + 1) // 2 + 1

def _total_index(x, y, cols, first):
    assert first == "up"
    if y % 2 == 0:
        return x + y // 2 * (_n_vert_cols(0, cols, first) + _n_vert_cols(1, cols, first))
    return x + y // 2 * _n_vert_cols(1, cols, first) + _ceildiv(y, 2) * _n_vert_cols(0, cols, first)

def _nearest_neighbours(x, y, first):
    assert first == "up"
    common = ((x, y-1),
              (x-1, y),
              (x, y+1),
              (x+1, y))
    if y % 2 == 0:
        return common + ((x+1, y+1), (x+1, y-1))
    return common + ((x-1, y-1), (x-1, y+1))

def _position(x, y, first, spacing):
    assert first == "up"
    if y % 2 == 0:
        return np.array(((x + 0.5) * spacing, np.sqrt(3/4)*y*spacing, 0))
    return np.array((x*spacing, np.sqrt(3/4)*y*spacing, 0))
