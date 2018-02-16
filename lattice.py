"""
Representation of a lattice.
"""

from copy import deepcopy
import numpy as np

def angle(p1, p2):
    "Compute angle between p1 and p2; result is in (-pi, pi]"
    return (np.arctan2(p2[1], p2[0]) - np.arctan2(p1[1], p1[0]) + np.pi) % (2*np.pi) - np.pi

def select_innermost(pos, centre):
    "Select site clostest to centre."
    return np.argmin([np.linalg.norm(centre-p) for p in pos])

def select_anticlockwise(pos, curpos):
    "Select next site in anti-clockwise order."
    return np.argmax([angle(curpos, p) for p in pos])

def select_next(method, pos, cursite, centre):
    "Select next site from list of positions."

    if method == "innermost":
        return select_innermost(pos, centre)
    if method == "anticlockwise":
        return select_anticlockwise(pos, cursite.pos)
    raise KeyError("Unknown traveral method: "+method+". Supported methods are:\
 innermost, anticlockwise")

def shifted_lattice(lat, vec):
    "Retrun a copy of a lattice shifted by a given vector."
    shifted = deepcopy(lat)
    for site in shifted:
        site.pos += vec
    return shifted


class Site:
    """
    A single lattice site.

    Attributes:
        idx - Index of the site.
        pos - (np.array) 2D or 3D position of the site.
        neighbours - List of indices of neighbouring sites.
        hopping - List of hopping strengths to sites listed in neighbours.
        Additional attributes can be stored using the member dict _attributes.
    """

    def __init__(self, idx, pos, neighbours=[], hopping=[], **attrs):
        self.idx = idx
        self.pos = pos
        self.neighbours = neighbours
        self.hopping = hopping
        self._attributes = attrs

    def has_attribute(self, attr):
        "Return True if this Site has an attribute called attr; False otherwise."
        return attr in self._attributes.keys()

    def __getitem__(self, attr):
        "Return an attribute."
        return self._attributes[attr]

    def __setitem__(self, attr, value):
        "Set an attribute."
        self._attributes[attr] = value

    def attributes(self):
        "Return an iterator over key, value pairs of attributes."
        return self._attributes.items()

    def relabel(self, labels):
        """
        Reassign site indices from a list of labels.
        Index i is replaced by index labels[i].
        """
        self.idx = labels[self.idx]
        self.neighbours = [labels[neigh] for neigh in self.neighbours]

class Lattice:
    "Represents a lattice."

    def __init__(self, name="", comment="", nt=0):
        self.name = name
        self.comment = comment
        self.nt = nt
        self.sites = []

    def __len__(self):
        "Number of sites."
        return len(self.sites)

    def __iter__(self):
        "Iterates over sites."
        return iter(self.sites)

    def __setitem__(self, attr, value):
        "Set an attribute for all sites in the lattice."
        for site in self.sites:
            site[attr] = value

    def at(self, pos, PREC=1e-13):
        "Return the site at a given position."
        for site in self.sites:
            if np.max(np.abs(pos-site.pos)) < PREC:
                return site
        return None

    def centre(self):
        "Get centre position of lattice."
        return sum(site.pos for site in self.sites)/len(self.sites)

    def check_consistency(self):
        "Check whether all site indices are in bounds."

        idxs = []
        # check site.idx
        for site in self.sites:
            if site.idx >= len(self.sites) or site.idx < 0:
                print("Site index out of range: {}".format(site.idx))
                return False
            idxs.append(site.idx)

        # check site.neighbours
        for site in self.sites:
            for neigh in site.neighbours:
                if neigh not in idxs:
                    print("Neighbour index not in registered site indexed: {}".format(neigh))
                    return False

        return True

    def label_graph(self, method):
        """
        Generate a list of labels (indices) for lattice sites. The labels spiral
        outward from the centre.

        This function is designed on 2D lattices and might not work in other dimensions!
        """

        cen = self.centre()
        # distance to centre for each site
        radii = [np.linalg.norm(cen-site.pos) for site in self.sites]
        # remember which sites were already visited
        visited = [False]*len(self.sites)

        # new labels
        labels = np.empty(len(self.sites), dtype=int)
        # start at position closest to centre
        cur = np.argmin(radii)
        labels[cur] = 0
        visited[cur] = True
        for i in range(1, len(self.sites)):
            try:
                # get all indices and positions of neighbours of cur that have not been visited
                idx, pos = zip(*[(j, self.sites[n].pos)
                                 for j, n in enumerate(self.sites[cur].neighbours)
                                 if not visited[n]])
            except ValueError:
                raise RuntimeError("Graph search got stuck. No neighbours to continue\
 after labeling {} sites".format(i+1))

            # go to next site
            cur = self.sites[cur].neighbours[idx[select_next(method, pos,
                                                             self.sites[cur], cen)]]
            labels[cur] = i
            visited[cur] = True
        return labels

    def relabel(self, labels):
        """
        Replace site indices by indices given in list labels.
        Index i is replaced by index labels[i].
        Also sorts the sites according to new indices.
        """
        for site in self.sites:
            site.relabel(labels)
        self.sites.sort(key=lambda site: site.idx)
