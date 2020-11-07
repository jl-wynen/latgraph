"""
Generator for pyramid 2D-tilings.
"""

r"""
IMPLEMENTATION NOTE
===================

The following shows how sites are indexed.
The numbers show the total indices while rows and columns are indexed
top-down and left-right, respectively.

       <-columns->
       0--- 1--- 2   base-row
       |\  /|\  /|
    ^  |  3 |  4 |   tip-row
    |  |/  \|/  \|
  rows 5--- 6--- 7   base-row
    |  |\  /|\  /|
    v  |  8 |  9 |   tip-row
       |/  \|/  \|
      10---11---12   base-row


nx and ny are the numbers of pyramids, not rows or columns.
The algorithm iterates through rows first, making for 2*ny+1 rows.
Then, columns are iterated, their number depends on the row and is nx+1 for
base-rows and nx for tip-rows.

Borders are identified based on the orientation in the figure above.
"t" = top
"b" = bottom
"l" = left
"r" = right
"""


import argparse
from enum import Enum, auto
import textwrap

import numpy as np

from lattice import Site, Lattice

def _define_parser():
    "Define command line argument parser for PyramidGen."

    parser = argparse.ArgumentParser(prog="latgraph --generate pyramid",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""\
                Generate a pyramid tiling in 2-D.

                Repeats an equilateral right square pyramid in two dimensions.
                Does not tile in the third dimension.

                Boundaries are open.
                                     """))
    parser.add_argument("nx", type=int, help="Number of pyramids in x direction.")
    parser.add_argument("ny", type=int, help="Number of pyramids in y direction.")
    parser.add_argument("--spacing", type=float, default=1.0, help="Lattice spacing.")
    parser.add_argument("--name", default="", help="Name for the lattice")
    parser.add_argument("--comment", default="", help="Comment on the lattice")
    return parser

def run(in_args):
    "Run PyramidGen from command line arguments."
    args = _define_parser().parse_args(in_args)
    return make_pyramid_tiling(args.nx, args.ny, args.spacing,
                               args.name, args.comment)

def make_pyramid_tiling(nx, ny, spacing, name, comment):
    """
    Create a Lattice by tiling pyramids.
    """

    lat = Lattice(name, comment)

    nsites = 2*nx*ny + nx + ny + 1
    i = 0
    base_row = True
    for row in range(2*ny + 1):
        for col in range(nx + (1 if base_row else 0)):
            pos = _site_position(row, col, base_row, spacing)
            neighbours = _neighbour_indices_of(i, base_row, nx,
                                               _borders(row, col, nx, ny, base_row))
            hoppings = [1]*len(neighbours)
            lat.sites.append(Site(i, pos, neighbours, hoppings))

            i += 1
        base_row = not base_row

    assert lat.check_consistency()
    return lat

def _borders(row, col, nx, ny, base_row):
    at_borders = ""
    if row == 0:
        at_borders += "t"
    elif row == 2*ny:
        at_borders += "b"
    if base_row:
        if col == 0:
            at_borders += "l"
        elif col == nx:
            at_borders += "r"
    return at_borders

def _neighbour_indices_of(i, base_row, nx, at_borders):
    def crosses_border(direction):
        return any(d in at_borders for d in direction)

    if base_row:
        index_shifts = (
            ("t", -2*nx - 1),
            ("tl", -nx - 1),
            ("l", -1),
            ("bl", nx),
            ("b", 2*nx + 1),
            ("br", nx + 1),
            ("r", 1),
            ("tr", -nx),
            )
    else:
        index_shifts = (
           ("tl", -nx - 1),
           ("bl", nx),
           ("br", nx + 1),
           ("tr", -nx),
        )

    return tuple(i+shift for direction, shift in index_shifts
                 if not crosses_border(direction))

def _site_position(row, col, base_row, spacing):
    if base_row:
        return np.array((col*spacing, row//2*spacing, 0))
    return np.array(((col + 0.5)*spacing, (row // 2 + 0.5)*spacing, spacing/np.sqrt(2)))
