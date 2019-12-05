"""
Generator for carbon nano tubes.
"""

import textwrap
import argparse
import numpy as np
import lattice

def _define_parser():
    "Define command line argument parser for AGNRgen."

    parser = argparse.ArgumentParser(prog="latgraph --generate agnr",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""\
                Generate an armchair carbon nanoribbon of given number of dimer lines (width) 
                and number of hexagons along the longitudinal direction (length).

                Currently only periodic boundary conditions (in the longitudinal direction) are supported. 
                                     """))
    parser.add_argument("dimer", type=int, help="Number of dimer lines that go in the longitudinal direction. Determines the width of the ribbon.")
    parser.add_argument("n_hex", type=int, help="Number of hexagons along the longitudinal direction. Determines the length of the ribbon.")
    parser.add_argument("--spacing", type=float, default="1", help="Lattice spacing")
    parser.add_argument("--name", default="", help="Name for the lattice")
    parser.add_argument("--comment", default="", help="Comment on the lattice")
    return parser


def run(in_args):
    "Run AGNRgen from command line arguments."

    parser = _define_parser()
    args = parser.parse_args(in_args)
	
    gen = AGNRgen(args.dimer, args.n_hex, args.spacing)

    lat = gen.make_agnr()
    #Use labels so that subplattice A has even and sublattice B odd indices
    lat.relabel(gen.bipartite_labels()) 
    
    if args.name:
        lat.name = args.name
    if args.comment:
        lat.comment = args.comment

    return lat

class AGNRgen:
    """
    Class to make armchair graphene nanoribbons (AGNRs), with variable width and length. 
    
    The width and length are determined via the following input arguments:
    dimer: number of dimer lines that go in the horizontal/longitudinal direction (determines width)
    n_hex: number of hexagons in the horizontal direction (determines length)
    spacing: lattice spacing
    
    Currently only periodic boundary conditions (along the horizontal/longitudinal direction) are implemented. 
    
    """
    def __init__(self, dimer, n_hex, spacing):
        self.dimer = dimer
        self.length = n_hex*2 #number of sites on one dimer line (in horizontal direction)
        self.spacing = spacing
    
    #Returns labels where one sublattice has even, the other has odd indices!
    def bipartite_labels(self): 
        N=self.dimer*self.length*2 #Is even!
        labels=[None]*N
        
        for idx in range(0, N):
            l, dim = self._get_l_dim(idx) 
            if dim%2==0:                #for even dimer lines nothing changes
                labels[idx]=idx
            else:                       #for odd dimer lines neighbours have to be flipped
                if l%2==0:
                    labels[idx]=idx+1
                else:
                    labels[idx]=idx-1
        return labels
          
    def unit_vectors(self):
        return np.array([3/2, np.sqrt(3)/2])*self.spacing, np.array([3/2, -np.sqrt(3)/2])*self.spacing
    
    def neighbour_vectors(self):
        #vectors that go from one site to all neighbours. Sign is sensitive on the "evenness" of the site.
        a1, a2 = self.unit_vectors()
        shift0 = 1/3 * (a1 + a2)
        shift1 = shift0 - a1
        shift2 = shift0 - a2
        return shift0, shift1, shift2
    
    #WORKS ONLY FOR INITIAL ORDERING!
    def _get_l_dim(self, idx):
        #l: place of the site in a dimer line, starting from the left most site as 0, can go from 0 to length
        #dim: dimmer line the site is one. Can go from 0 to dimer.
        dim = idx//self.length
        l=idx-self.length*dim
        return l, dim

    #WORKS ONLY FOR INITIAL ORDERING!
    def _get_neighbours(self, idx): 
        #returns indices of the neighbours as a list of a single site with given index.
        neigh=[]
        l, dim = self._get_l_dim(idx)
        
        #distiguish between the two sublattices
        if (dim%2+l)%2==0: #site is "even", meaning it has an upper right neighbour 
            #lower right neighbor:
            if dim>0:
                neigh.append(idx-self.length)
            #left neighbor:
            if l>0:
                neigh.append(idx-1)
            else:
                neigh.append(idx+self.length-1)
            #upper right neigbour:
            if dim<self.dimer-1:
                neigh.append(idx+self.length)
        else: #site is "odd", meaning it has an upper left neighbour 
            #lower left neighbor:
            if dim>0:
                neigh.append(idx-self.length)
            #right neighbor:
            if l<self.length-1:
                neigh.append(idx+1)
            else:
                neigh.append(idx-self.length+1)
            #upper left neighbor:
            if dim<self.dimer-1:
                neigh.append(idx+self.length)
        return neigh
    
    def make_agnr(self):
        #returns an AGNR lattice with periodic boundary conditions and all hopping strengths set to 1
        #and with indices ordered by starting in the bottom left and then going along the dimer lines
        #to the right. When the boundary start go to the left of the above dimer line.
        
        
        #geometry of lattice to determine the position of the sites
        unit_vectors = self.unit_vectors()
        shift0, shift1, shift2 = self.neighbour_vectors()
        distance_dim=np.abs(shift2[1])
        
        lat=lattice.Lattice()
        
        #Generate sites (with positions)
        for dim in range(0, self.dimer):
            xpos=(dim%2)*(-shift1[0])     #odd dimer lines start indented by 0.5*spacing, even start at 0
            for l in range(0, self.length):
                idx=dim*self.length+l
                ypos=distance_dim*dim
                lat.sites.append(lattice.Site(idx, np.array([xpos, ypos])))
                if (l+dim%2)%2==0:
                    xpos += (-shift1[0]+unit_vectors[1][0])
                else:
                    xpos += self.spacing
        
        #Fill in neighbours for each site with periodic boundary conditions and set hopping to 1
        for site in lat.sites:
            neigh=self._get_neighbours(site.idx)
            site.neighbours=neigh
            site.hopping=[1]*len(neigh)
        return lat
        