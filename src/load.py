
from __future__ import annotations

import numpy as np

from src.matrices import elementary_mass_matrix


def load_vector(vtx: np.ndarray, tri: np.ndarray, f_nodes: np.ndarray) -> np.ndarray:
    n = len(vtx)
    F = np.zeros(n)

    for t in tri:
        coords = vtx[t]
        Mloc = elementary_mass_matrix(coords)
        Floc = Mloc @ f_nodes[t]

        for i_local, i_global in enumerate(t):
            F[int(i_global)] += Floc[i_local]

    return F
