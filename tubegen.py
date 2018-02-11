"""
See Physical Properties Of Carbon Nanotubes
     -  Dresselhaus G,Dresselhaus Mildred S,Saito
for definitions of vectors
"""

from copy import deepcopy
from math import gcd

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import lattice

class TubeGen:
    def __init__(self, chirality, spacing):
        self.chirality = chirality
        self.spacing = spacing


    def n_hex_ucell(self):
        n, m = self.chirality
        return 2*(m*m + n*n + m*n)//gcd(2*m+n, 2*n+m)

    def n_atoms_ucell(self):
        return 2*self.n_hex_ucell()

    def circumference(self):
        return np.linalg.norm(self.chiral_vector())

    def diameter(self):
        return self.circumference()/np.pi

    def unit_vectors(self):
        return np.array([3/2, np.sqrt(3)/2])*self.spacing, \
            np.array([3/2, -np.sqrt(3)/2])*self.spacing

    def chiral_vector(self):
        a1, a2 = self.unit_vectors()
        return self.chirality[0]*a1 + self.chirality[1]*a2

    def translation_vector_lat_basis(self):
        n, m = self.chirality
        dR = gcd(2*m+n, 2*n+m)
        return np.array(((2*m+n)/dR, -(2*n+m)/dR))

    def translation_vector(self):
        a1, a2 = self.unit_vectors()
        t1, t2 = self.translation_vector_lat_basis()
        return t1*a1 + t2*a2
    
    def symmetry_vector(self):
        a1, a2 = self.unit_vectors()
        t1, t2 = self.translation_vector_lat_basis()

        # require t1*q - t2*p = 1
        # and find smallest p that fulfills this equation for integer p, q
        for p in range(1, self.n_atoms_ucell()+1):
            q = (1 + p*t2)/t1
            if np.abs(q%1.) < 1e-10:  # q close to integer
                return p*a1 + q*a2

        raise RuntimeError("Unable to find symmetry vector")

    def make_ucell(self):
        T = self.translation_vector()
        Ch = self.chiral_vector()
        R = self.symmetry_vector()
        t = np.linalg.norm(T)
        ch = np.linalg.norm(Ch)
        uCh = Ch / ch
        uT = T / t

        a1, a2 = self.unit_vectors()
        eo_shift = 1/3 * (a1 + a2)

        sites = []
        idx = 0
        for i in range(self.n_atoms_ucell()//2):
            # even
            psi = np.dot(i*R, uCh) % ch
            tau = np.dot(i*R, uT) % t
            sites.append(lattice.Site(idx, psi*uCh + tau*uT))
            idx += 1

            # odd
            psi = np.dot(i*R + eo_shift, uCh) % ch
            tau = np.dot(i*R + eo_shift, uT) % t
            sites.append(lattice.Site(idx, psi*uCh + tau*uT))
            idx += 1

        return sites

    def make_sites(self, n_ucells):
        ucell = self.make_ucell()
        T = self.translation_vector()
        sites = []
        for i in range(n_ucells):
            sites.extend(shift_lattice(ucell, T*i))
        return sites

    def __str__(self):
        a1, a2 = self.unit_vectors()
        Ch = self.chiral_vector()
        T = self.translation_vector()
        R = self.symmetry_vector()

        return """## TubeGen ##
Chirality:           ({ch1}, {ch2})
Hexes per unit cell: {nhpuc}
Atoms per unit cell: {napuc}
Circumference:       {circ}
Diameter:            {diam}
Unit vectors:        ({uv11}, {uv12}), ({uv21}, {uv22})
Chiral vector:       ({cv1}, {cv2})
Translation Vector:  ({tv1}, {tv2})
Symmetry Vector:     ({sv1}, {sv2})""".format(
    ch1=self.chirality[0], ch2=self.chirality[1],
    nhpuc=self.n_hex_ucell(),
    napuc=self.n_atoms_ucell(),
    circ=self.circumference(),
    diam=self.diameter(),
    uv11=a1[0], uv12=a1[1], uv21=a2[0], uv22=a2[1],
    cv1=Ch[0], cv2=Ch[1],
    tv1=T[0], tv2=T[1],
    sv1=R[0], sv2=R[1]
)

def shift_lattice(sites, vec):
    shifted = deepcopy(sites)
    for site in shifted:
        site.pos += vec
    return shifted

def rotate(vec, angle):
    r = np.linalg.norm(vec)
    phi = np.arctan2(vec[1], vec[0])
    return r*np.array((np.cos(phi-angle), np.sin(phi-angle)))

def rotate_lattice(sites, angle):
    rotated = deepcopy(sites)
    for site in rotated:
        site.pos = rotate(site.pos, angle)
    return rotated

def main():
    gen = TubeGen((4, 2), 1)
    sites = gen.make_sites(2)
    Ch = gen.chiral_vector()
    # sites = rotate_lattice(sites, np.arctan2(Ch[1], Ch[0]))

    print(gen)


    even = []
    odd = []
    for i, site in enumerate(sites):
        if i % 2 == 0:
            even.append(site.pos)
        else:
            odd.append(site.pos)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.scatter(*zip(*even), label="even")
    ax.scatter(*zip(*odd), label="odd")

    for site in sites:
        ax.text(*site.pos, str(site.idx))
    Ch = gen.chiral_vector()
    T = gen.translation_vector()
    R = gen.symmetry_vector()
    ax.plot((0, Ch[0]), (0, Ch[1]), label="Ch")
    ax.plot((0, T[0]), (0, T[1]), label="T")
    ax.plot((0, R[0]), (0, R[1]), label="R")

    ax.legend()
    fig.tight_layout()

    plt.show()



if __name__ == "__main__":
    main()
