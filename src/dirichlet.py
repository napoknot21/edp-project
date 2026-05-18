
from __future__ import annotations

import numpy as np
import scipy.sparse.linalg as spla


def boundary_vertices (bnd: np.ndarray) -> np.ndarray :
    """
    
    """
    return np.unique (bnd.reshape(-1))


def solve_zero_dirichlet (A, rhs, bnd: np.ndarray, n_vertices: int) -> np.ndarray:
    """

    """
    bd = boundary_vertices(bnd)
    all_dofs = np.arange(n_vertices)
    free = np.setdiff1d(all_dofs, bd)

    A_free = A[free, :][:, free].tocsr()
    rhs_free = rhs[free]

    u = np.zeros(n_vertices)
    u[free] = spla.spsolve(A_free, rhs_free)

    return u


def solve_dirichlet_values(A, rhs, values: np.ndarray, constrained: np.ndarray) -> np.ndarray :
    """
    """
    n = A.shape[0]
    all_dofs = np.arange(n)
    free = np.setdiff1d(all_dofs, constrained)

    u = np.zeros(n)
    u[constrained] = values[constrained]

    rhs_free = rhs[free] - A[free, :][:, constrained] @ u[constrained]
    A_free = A[free, :][:, free].tocsr()

    u[free] = spla.spsolve(A_free, rhs_free)
    return u
