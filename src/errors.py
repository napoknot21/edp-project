
from __future__ import annotations

import numpy as np


def l2_norm (u, M) -> float :
    """
    Compute the discrete L2 norm associated with the mass matrix.
    """
    return float(np.sqrt(u @ (M @ u)))


def h1_norm (u, M, K) -> float :
    """
    Compute the discrete H1 norm associated with mass and stiffness matrices.
    """
    return float(np.sqrt(u @ (M @ u) + u @ (K @ u)))


def relative_l2_error (uh, uref, M) -> float :
    """
    Compute the relative L2 error between a numerical and reference vector.
    """
    e = uh - uref
    return l2_norm(e, M) / l2_norm(uref, M)


def relative_h1_error (uh, uref, M, K) -> float :
    """
    Compute the relative H1 error between a numerical and reference vector.
    """
    e = uh - uref
    return h1_norm(e, M, K) / h1_norm(uref, M, K)


def mesh_size (vtx, tri) -> float :
    """
    Return the largest triangle edge length in the mesh.
    """
    h = 0.0

    for t in tri :
    
        pts = vtx[t]
    
        h = max(
    
            h,
            np.linalg.norm(pts[1] - pts[0]),
            np.linalg.norm(pts[2] - pts[1]),
            np.linalg.norm(pts[0] - pts[2]),
    
        )
    
    return float(h)
