
from __future__ import annotations

import numpy as np


def triangle_area(coords: np.ndarray) -> float:
    x0, y0 = coords[0]
    x1, y1 = coords[1]
    x2, y2 = coords[2]
    return 0.5 * abs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0))


def edge_length(coords: np.ndarray) -> float:
    return float(np.linalg.norm(coords[1] - coords[0]))


def p1_gradients(coords: np.ndarray) -> np.ndarray:
    x = coords[:, 0]
    y = coords[:, 1]

    two_area_signed = (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0])

    if abs(two_area_signed) < 1e-14:
        raise ValueError("Degenerate triangle")

    grads = np.zeros((3, 2))

    grads[0] = [y[1]-y[2], x[2]-x[1]]
    grads[1] = [y[2]-y[0], x[0]-x[2]]
    grads[2] = [y[0]-y[1], x[1]-x[0]]

    grads /= two_area_signed
    return grads
