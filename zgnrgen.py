"""
Generator for carbon nano ribbons with zig-zag edges (ZGNRs).
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
                Generate a zig-zag carbon nanoribbon of given number of zig-zag lines (width) 
                and number of unit cells along the longitudinal direction (length).

                Currently only periodic boundary conditions are supported (in the longitudinal direction). 
                                     """))
    parser.add_argument("n_zig", type=int, help="Number of dimer lines that go in the longitudinal direction. Determines the width of the ribbon.")
    parser.add_argument("n_uc", type=int, help="Number of hexagons along the longitudinal direction. Determines the length of the ribbon.")
    parser.add_argument("--spacing", type=float, default="1", help="Lattice spacing")
    parser.add_argument("--name", help="Name for the lattice. Default is N_ZIGzgnrN_UC")
    parser.add_argument("--comment", default="", help="Comment on the lattice")
    return parser


def run(in_args):
    "Run ZGNRgen from command line arguments."
    parser = _define_parser()
    args = parser.parse_args(in_args)
	
    gen = ZGNRgen(args.n_zig, args.n_uc, args.spacing)

    lat = gen.make_lattice()

    if args.name:
        lat.name = args.name
    else:
        lat.name = "{}zgnr{}".format(args.n_zig, args.n_uc)
    if args.comment:
        lat.comment = args.comment

    return lat

class ZGNRgen:
    """
    Class to make zig-zag graphene nanoribbons (ZGNRs), with variable width and length. 
    
    The width and length are determined via the following input arguments:
    n_zig: number of zig-zag lines (determines width)
    n_uc: number of unit cells in the longitudinal direction (determines length)
    spacing: lattice spacing
    
    Currently only periodic boundary conditions (along the longitudinal direction) are implemented. 
    Boundary conditions along the transverse direction are open, otherwise we would call it a tube.
    """
    def __init__(self, n_zig, n_uc, spacing):
        self.n_zig = n_zig #Number of zigzag lines 
        self.n_uc = n_uc #Number of hexagons in longitudinal direction
        self.spacing = spacing
    
   
    def get_coordinates(self, idx):
        """Returns:
            - zig: The index of the zig-zag line, counted from the bottom up)
            - l: The index in the zig-zag line counted as
            zig even:
                1       3
            0       2       4
            zig odd:
            1       3       
                0       2
            """
        zig=idx//(2*self.n_uc)
        l=idx%(self.n_uc*2)
        return zig, l
  
    def get_neighbours(self, idx):
        """"Return the indices of all neighbours in an array, applying periodic
        boundary conditions"""
        
        zig, l = self.get_coordinates(idx)
        neigh=[]
        
        lmax=self.n_uc*2-1
        
        if idx%2==0: #Site has two upper neighbors
            #lower neighbour:
            if zig>0:               #Site isn't on the bottom boundary
                neigh.append(idx-lmax)
            
            #upper left:
            if zig%2==0 and l==0:   #Periodic boundary
                neigh.append(idx+lmax)
            else:
                neigh.append(idx+(-1+2*(zig%2)))
            
            #upper right:
            if zig%2==1 and (l==lmax-1):   #Periodic boundary
                neigh.append(idx-lmax+2)
            else:
                neigh.append(idx+1+(zig%2)*2)
           
  
        else:   #Site has two lower neighbours   
             #lower left:
            if zig%2==1 and l==1: #Periodic boundary
                neigh.append(idx+lmax-2)
            else:
                neigh.append(idx-1-(zig%2)*2)
            
            #lower right:
            if zig%2==0 and l==lmax:
                neigh.append(idx-lmax)
            else:
                neigh.append(idx+1-2*(zig%2))
                
            #upper neighbor:
            if zig<(self.n_zig-1):       #Site isn't on the top boundary
                neigh.append(idx+lmax)
        return neigh
  
    
    def make_lattice(self):        
        #Generate sites (with positions, but no neighbours)
        lat=lattice.Lattice()
       
        for zig in range(0, self.n_zig):
            if zig%2==0:
                for l in range(self.n_uc*2):
                    xpos=l//2*np.sqrt(3) + np.sqrt(3)/2*(l%2)
                    ypos=zig*3/2+1/2*(l%2)
                    idx1 = zig*self.n_uc*2+l
                    lat.sites.append(lattice.Site(idx1, np.array([xpos*self.spacing, ypos*self.spacing])))
            else:
                for l in range(self.n_uc*2):
                    xpos=l//2*np.sqrt(3) + np.sqrt(3)/2*((l+1)%2)
                    ypos=zig*3/2+1/2*(l%2)
                    idx1 = zig*self.n_uc*2+l
                    lat.sites.append(lattice.Site(idx1, np.array([xpos*self.spacing, ypos*self.spacing])))
        
        #Fill in neighbours for each site with periodic boundary conditions and set hopping to 1
        for site in lat.sites:
            neigh=self.get_neighbours(site.idx)
            site.neighbours=neigh
            site.hopping=[1]*len(neigh)
        return lat