
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from src.geometry import triangle_area, edge_length, p1_gradients


def elementary_mass_matrix(coords: np.ndarray) -> np.ndarray:
    area = triangle_area(coords)
    return area / 12.0 * np.array([
        [2, 1, 1],
        [1, 2, 1],
        [1, 1, 2],
    ], dtype=float)


def elementary_edge_mass_matrix(coords: np.ndarray) -> np.ndarray:
    length = edge_length(coords)
    return length / 6.0 * np.array([
        [2, 1],
        [1, 2],
    ], dtype=float)


def elementary_stiffness_matrix(coords: np.ndarray) -> np.ndarray:
    area = triangle_area(coords)
    grads = p1_gradients(coords)
    return area * (grads @ grads.T)


def elementary_convection_matrix(coords: np.ndarray, b: np.ndarray) -> np.ndarray:
    area = triangle_area(coords)
    grads = p1_gradients(coords)

    local = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            local[i, j] = area / 3.0 * np.dot(b, grads[j])

    return local


def assemble_triangle_matrix(vtx: np.ndarray, tri: np.ndarray, elementary_function, *args) -> sp.csr_matrix:
    n = len(vtx)
    rows, cols, data = [], [], []

    for t in tri:
        coords = vtx[t]
        local = elementary_function(coords, *args)

        for i_local, i_global in enumerate(t):
            for j_local, j_global in enumerate(t):
                rows.append(int(i_global))
                cols.append(int(j_global))
                data.append(float(local[i_local, j_local]))

    return sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()


def mass_matrix(vtx: np.ndarray, tri: np.ndarray) -> sp.csr_matrix:
    return assemble_triangle_matrix(vtx, tri, elementary_mass_matrix)


def stiffness_matrix(vtx: np.ndarray, tri: np.ndarray) -> sp.csr_matrix:
    return assemble_triangle_matrix(vtx, tri, elementary_stiffness_matrix)


def convection_matrix(vtx: np.ndarray, tri: np.ndarray, b: np.ndarray) -> sp.csr_matrix:
    return assemble_triangle_matrix(vtx, tri, elementary_convection_matrix, b)


def boundary_mass_matrix(vtx: np.ndarray, bnd: np.ndarray) -> sp.csr_matrix:
    n = len(vtx)
    rows, cols, data = [], [], []

    for e in bnd:
        coords = vtx[e]
        local = elementary_edge_mass_matrix(coords)

        for i_local, i_global in enumerate(e):
            for j_local, j_global in enumerate(e):
                rows.append(int(i_global))
                cols.append(int(j_global))
                data.append(float(local[i_local, j_local]))

    return sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
