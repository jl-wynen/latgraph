"""

"""

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

def main():
    lat3 = fileio.read_w3d("c60_ipr.w3d")
    if not lat3.check_consistency():
        return
    lat2 = fileio.read_w2d("c60_ipr.w2d")
    if not lat2.check_consistency():
        return

    show_lattice(lat3, "Original", np.arange(len(lat3)))

    labels = lat2.label_graph()
    lat3.relabel(labels)
    show_lattice(lat3, "Relabelled", np.arange(len(lat3)))

    fileio.write_wnd("c60_ipr_rel.w3d", lat3)
    fileio.write_wnd("c60_ipr_rel.w2d", lat2)

    plt.show()

if __name__ == "__main__":
    main()
