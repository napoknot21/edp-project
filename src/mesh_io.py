
from __future__ import annotations

import numpy as np
import meshio


def load_mesh(path: str):
    mesh = meshio.read(path)
    vtx = np.asarray(mesh.points[:, :2], dtype=float)

    tri = np.asarray(mesh.cells_dict.get("triangle", []), dtype=int)

    if "line" in mesh.cells_dict:
        bnd = np.asarray(mesh.cells_dict["line"], dtype=int)
    else:
        bnd = extract_boundary_edges(tri)

    return vtx, tri, bnd


def extract_boundary_edges(tri: np.ndarray) -> np.ndarray:
    edge_count = {}

    for t in tri:
        edges = [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]
        for a, b in edges:
            e = tuple(sorted((int(a), int(b))))
            edge_count[e] = edge_count.get(e, 0) + 1

    bnd = [e for e, count in edge_count.items() if count == 1]
    return np.asarray(bnd, dtype=int)
