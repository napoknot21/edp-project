
from __future__ import annotations

import numpy as np


def source_values(vtx: np.ndarray, p: int, q: int, b: np.ndarray) -> np.ndarray:
    x = vtx[:, 0]
    y = vtx[:, 1]

    return (
        np.exp(0.5 * (b[0] * x + b[1] * y))
        * np.sin(p * np.pi * x)
        * np.sin(q * np.pi * y)
    )


def alpha_value(p: int, q: int, b: np.ndarray, c: float) -> float:
    denominator = np.pi**2 * (p**2 + q**2) + c + 0.25 * np.dot(b, b)
    return float(1.0 / denominator)


def exact_values(vtx: np.ndarray, p: int, q: int, b: np.ndarray, c: float) -> np.ndarray:
    return alpha_value(p, q, b, c) * source_values(vtx, p, q, b)
