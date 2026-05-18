
from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


def plot_mesh (vtx, tri, output) :
    """
    Save a PNG figure of the triangular mesh.
    """
    output = Path(output)
    
    plt.figure(figsize=(6, 6))
    plt.triplot(vtx[:, 0], vtx[:, 1], tri, linewidth=0.5)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Triangular mesh")
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(output, dpi=200)
    plt.close()


def plot_solution (vtx, tri, values, output, title) :
    """
    Save a filled contour plot of nodal values on a triangular mesh.
    """
    output = Path(output)
    
    triang = mtri.Triangulation(vtx[:, 0], vtx[:, 1], tri)

    plt.figure(figsize=(7, 6))
    c = plt.tricontourf(triang, values, levels=40)
    
    plt.colorbar(c, label="value")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(output, dpi=200)
    
    plt.close()


def plot_convergence(h, e_l2, e_h1, output):
    """
    Save a log-log convergence plot for relative L2 and H1 errors.
    """
    output = Path(output)

    order = np.argsort(h)
    h = h[order]
    e_l2 = e_l2[order]
    e_h1 = e_h1[order]

    plt.figure(figsize=(7, 5))
    plt.loglog(h, e_l2, "o-", label="relative L2 error")
    plt.loglog(h, e_h1, "s-", label="relative H1 error")

    if len(h) >= 2 :

        plt.loglog(h, e_l2[-1] * (h / h[-1])**2, "--", label="slope 2 reference")
        plt.loglog(h, e_h1[-1] * (h / h[-1]), "--", label="slope 1 reference")

    plt.xlabel("mesh size h")
    plt.ylabel("relative error")
    plt.title("Convergence")
    plt.legend()
    plt.grid(True, which="both")
    plt.tight_layout()
    plt.savefig(output, dpi=200)
    
    plt.close()
