"""
Representation of a lattice.
"""

import numpy as np

class Site:
    """
    A single lattice site.

    Attributes:
        idx - Index of the site.
        pos - (np.array) 2D or 3D position of the site.
        neighbours - List of indices of neighbouring sites.
        hopping - List of hopping strengths to sites listed in neighbours.
    """

    def __init__(self, idx, pos, neighbours, hopping):
        self.idx = idx
        self.pos = pos
        self.neighbours = neighbours
        self.hopping = hopping

    def relabel(self, labels):
        """
        Reassign site indices from a list of labels.
        Index i is replaced by index labels[i].
        """
        self.idx = labels[self.idx]
        self.neighbours = [labels[neigh] for neigh in self.neighbours]

class Lattice:
    "Represents a lattice."

    def __init__(self):
        self.sites = []

    def __len__(self):
        "Number of sites."
        return len(self.sites)

    def __iter__(self):
        "Iterates over sites."
        return iter(self.sites)

    def centre(self):
        "Get centre position of lattice."
        return sum(site.pos for site in self.sites)

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

    def label_graph(self):
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
                # get all indices and radii of neighbours of cur that have not been visited
                idx, rad = zip(*[(j, radii[n])
                                 for j, n in enumerate(self.sites[cur].neighbours)
                                 if not visited[n]])
            except ValueError:
                raise RuntimeError("Graph search got stuck. No neighbours to continue\
 after labeling {} sites".format(i+1))

            # move to site closest to centre
            cur = self.sites[cur].neighbours[idx[np.argmin(rad)]]
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
