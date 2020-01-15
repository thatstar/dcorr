import numpy as np
from ase import Atoms
from pyasa.voronoi import VoronoiAnalysis


def voronoi_indices(snapshot):
    atoms = Atoms(cell=snapshot["box"],
                  positions=snapshot["positions"].copy(),
                  pbc=True)
    atoms.wrap()
    voro = VoronoiAnalysis(atoms)
    d = voro.get_voronoi_statistics()
    inds = np.empty((len(atoms), 4), dtype=int)
    for i in d.keys():
        for j in d[i]:
            inds[j, 0] = i[2]
            inds[j, 1] = i[3]
            inds[j, 2] = i[4]
            inds[j, 3] = i[5]

    return inds
