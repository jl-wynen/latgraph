######## Generate a Kagome Lattice #########

import sys
import numpy as np
import matplotlib.pyplot as plt
import yaml
from lattice import Lattice, Site
import argparse
import textwrap

# kagome lattice unit cell has three particles let them be  at [a(0,0), b(1,0), c(0.5,0.866)]
# basis vectors [(1,0), (0.5,0.866)]
# creates only 2-d lattice


def _define_parser():
    "Define command line argument parser for KagomeLatticeGen."
    parser = argparse.ArgumentParser(prog="latgraph --generate kagome",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""\
                Generate a kagome lattice.
                                     """))
    parser.add_argument("cols", type=int, help="Number of unit cell per row.")
    parser.add_argument("rows", type=int, help="Number of stacked rows.")
    parser.add_argument("boundary_condition", type=int, choices=(1, 0),
                        help="Boundary condition --periodic=1,open=0")
    parser.add_argument("--spacing", type=float, default=1.0, help="Lattice spacing.")
    parser.add_argument("--name", default="", help="Name for the lattice")
    parser.add_argument("--comment", default="", help="Comment on the lattice")
    return parser


def run(in_args):
    "Run Kagome lattice from command line arguments."

    args = _define_parser().parse_args(in_args)
    if args.boundary_condition == 1:
        print("Lattice with periodic boundary conditions")
    else:
        print("Lattice with out periodic boundary conditions")

    return make_kagome(args.cols, args.rows, args.boundary_condition,args.spacing,
                       args.name, args.comment)


def position(Lx, Ly, trans):

    a = []  # position of lattice a
    b = []  # position of lattice b
    c = []  # position of lattice c
    sites = []  # position of all lattice sites in the order a,b,c
    s = 0
    size = Lx*Ly

    for i in range(size):
        if i % Lx == 0 and i != 0:
            s = s+2
        a.append([0 + (i % Lx+s)*trans[0], 0 + i % Ly*trans[1]])
        sites.append(a[-1])
        b.append([1 + (i % Lx+s)*trans[0], 0 + i % Ly*trans[1]])
        sites.append(b[-1])
        c.append([0.5 + (i % Lx+s)*trans[0], 0.5*np.sqrt(3) + i % Ly*trans[1]])
        sites.append(c[-1])
    return a, b, c, sites


def nearest_neighbours(size, a, b, c):
    # arrays to get nearest neighbour connections
    ne_a = np.array(
        [[1.0, 0.0], [-1.0, 0.0], [0.5, np.sqrt(3)*0.5], [-0.5, -np.sqrt(3)*0.5]])
    ne_b = np.array(
        [[1.0, 0.0], [-1.0, 0.0], [0.5, -np.sqrt(3)*0.5], [-0.5, np.sqrt(3)*0.5]])
    ne_c = np.array([[0.5, np.sqrt(3)*0.5], [-0.5, np.sqrt(3)*0.5], [0.5, -np.sqrt(3)*0.5],
                     [-0.5, -np.sqrt(3)*0.5]])

    # nearest neighbour positions
    n_pos_a = [[i + ne_a] for i in a]
    n_pos_b = [[i + ne_b] for i in b]
    n_pos_c = [[i + ne_c] for i in c]

    # for each element in sites sites_pos gives its nearest neighbours
    site_pos = []
    for i in range(size):
        site_pos.append(n_pos_a[i])
        site_pos.append(n_pos_b[i])
        site_pos.append(n_pos_c[i])

    return site_pos


def _finding_index(boundary,sites):
    ind = []
    for i in boundary:
        # finding the index of the boundary sites in sites[]
        us = list(np.where(np.round(sites, 5) == np.round(i, 5)))[0]
        # print(us)
        res = [idx for idx, val in enumerate(us) if val in us[:idx]]
        # print(res)
        if res != []:
            ind.append(us[res[0]])
    return ind


def adjacencyList(Lx, Ly, site_pos, sites, boundary_ind_a, boundary_ind_b, boundary_ind_c, periodic):
    neighbour_index =[]
    adjacency_list = []
    v1 = v2 = v3 = v4 = v5 = v6 = 0
    for i in range(len(site_pos)):
        ind = []
        for j in range(len(site_pos[i][0])):
            # finding the index of the neighbour sites in sites[]
            us = list(np.where(np.round(sites, 5) ==
                               np.round(site_pos[i][0][j], 5)))[0]
            # print(us)
            res = [idx for idx, val in enumerate(us) if val in us[:idx]]
            # print(res)
            if res != []:
                ind.append(us[res[0]])

        if periodic == 1:
            if len(ind) != 4:
                if i in boundary_ind_b:
                    if i < (Lx*3*(Ly-1)+1):
                        ind.append(boundary_ind_c[Ly+v1])
                        v1 = v1+1
                    else:
                        ind.append(boundary_ind_c[Ly+v1])
                        ind.append(boundary_ind_a[Lx+v2])
                        v1 = -Ly + v2
                        v2 = v2+1
                elif i in boundary_ind_a:
                    if i <= (Ly*3-1):
                        ind.append(boundary_ind_b[Ly+v3])
                        v3 = v3+1
                    if (i > (Ly*3-1) or i == 0):
                        ind.append(boundary_ind_c[Lx+v4])
                        v4 = v4+1
                else:
                    if i < (Ly*3-1):
                        ind.append(boundary_ind_b[Ly+1+v5])
                        v5 = v5+1
                    else:
                        ind.append(boundary_ind_a[v6])
                        ind.append(boundary_ind_b[v6])
                        v6 = v6+1
            ind.sort()
            neighbour_index.append(ind)
        else:
            ind.sort()
            neighbour_index.append(ind)
            # if len(ind) !=4:
            #     print('bad')
            # else:
            #     print(i,ind)

        for k in ind:
            if [i, k] and [k, i] not in adjacency_list:  # removing symmetery
                adjacency_list.append([i, k])

    adjacency_list = np.array(adjacency_list)
    return neighbour_index, adjacency_list


def make_kagome(Lx, Ly, periodic,spacing, name, comment):

    lat = Lattice(name, comment)
    size = Lx*Ly
    num_sites = size*3

    # translation vector
    trans =list(np.array([1.0, np.sqrt(3)])*spacing)
    a, b, c, sites = position(Lx, Ly, trans)
    site_pos = nearest_neighbours(size, a, b, c)
    # boundary sites
    boundary_sites_a = sites[0::Lx*3]+sites[0:Ly*3:3]
    boundary_sites_b = sites[1::Lx*3]+sites[-(Ly*3-1)::3]
    boundary_sites_c = sites[2:(Ly*3):3]+sites[(Ly*3)-1::(Lx*3)]

    boundary_ind_a = _finding_index(boundary_sites_a,sites)
    boundary_ind_b = _finding_index(boundary_sites_b,sites)
    boundary_ind_c = _finding_index(boundary_sites_c,sites)

    neighbour_index, adjacency = adjacencyList(Lx, Ly,
                              site_pos, sites, boundary_ind_a, boundary_ind_b, boundary_ind_c, periodic)
   
    #hopping_matrix = np.ones(len(adjacency))
    positions = np.zeros((num_sites, 3))
    positions[:, :-1] = sites  # making 3-dim
    for i in range(len(sites)):
        hoppings = [1]*len(neighbour_index[i]) # hopping matrix of equal amplitude
        lat.sites.append(Site(i,positions[i],neighbour_index[i],hoppings))

    return lat
