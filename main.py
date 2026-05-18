
from __future__ import annotations

import argparse
import numpy as np

from pathlib import Path

from src.mesh_io import load_mesh
from src.matrices import mass_matrix, stiffness_matrix, convection_matrix
from src.load import load_vector
from src.dirichlet import solve_zero_dirichlet
from src.exact import source_values, exact_values, alpha_value
from src.errors import relative_l2_error, relative_h1_error, mesh_size
from src.plotting import plot_mesh, plot_solution, plot_convergence


def solve_one (mesh_path, p, q, bx, by, c, out, plot=False) :
    
    vtx, tri, bnd = load_mesh(mesh_path)
    b = np.array([bx, by], dtype=float)

    print(f"\n[*] Mesh: {mesh_path}")
    print(f"[*] Number of vertices: {len(vtx)}")
    print(f"[*] Number of triangles: {len(tri)}")
    print(f"[*] Number of boundary edges: {len(bnd)}")
    print(f"[+] h = {mesh_size(vtx, tri):.6e}")

    M = mass_matrix(vtx, tri)
    K = stiffness_matrix(vtx, tri)
    C = convection_matrix(vtx, tri, b)

    A = K + C + c * M

    f_nodes = source_values(vtx, p, q, b)
    F = load_vector(vtx, tri, f_nodes)

    uh = solve_zero_dirichlet(A, F, bnd, len(vtx))
    uex = exact_values(vtx, p, q, b, c)

    e_l2 = relative_l2_error(uh, uex, M)
    e_h1 = relative_h1_error(uh, uex, M, K)

    print(f"\n[*] alpha = {alpha_value(p, q, b, c):.10e}")
    print(f"[*] relative L2 error = {e_l2:.6e}")
    print(f"[*] relative H1 error = {e_h1:.6e}")

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)

    mesh_name = Path(mesh_path).stem

    if plot :

        plot_mesh(vtx, tri, out / f"mesh_{mesh_name}.png")
        plot_solution(vtx, tri, uh, out / f"uh_{mesh_name}.png", "Numerical solution $u_h$")
        plot_solution(vtx, tri, uex, out / f"uex_nodes_{mesh_name}.png", "Exact solution at nodes")
        
        print(f"\n[+] Figures saved in {out}")

    return mesh_size(vtx, tri), e_l2, e_h1


def run_convergence(args) :
    """
    
    """
    hs, e_l2s, e_h1s = [], [], []
    
    for mesh_path in args.meshes :

        h, e_l2, e_h1 = solve_one(
            mesh_path, args.p, args.q, args.bx, args.by, args.c, args.out, plot=False
        )
        hs.append(h)
        
        e_l2s.append(e_l2)
        e_h1s.append(e_h1)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    
    tag = f"p{args.p}_q{args.q}_bx{args.bx:g}_by{args.by:g}_c{args.c:g}"

    plot_convergence(np.array(hs), np.array(e_l2s), np.array(e_h1s), out / f"convergence_{tag}.png")
    
    print(f"[*] Convergence plot saved in {out / f'convergence_{tag}.png'}")


def parser() :

    p = argparse.ArgumentParser()
    
    p.add_argument("--mode", choices=["single", "convergence"], default="single")
    p.add_argument("--mesh", type=str)
    p.add_argument("--meshes", nargs="*")
    p.add_argument("--p", type=int, default=1)
    p.add_argument("--q", type=int, default=1)
    p.add_argument("--bx", type=float, default=1.0)
    p.add_argument("--by", type=float, default=1.0)
    p.add_argument("--c", type=float, default=1.0)
    p.add_argument("--out", type=str, default="outputs")
    p.add_argument("--plot", action="store_true")
    
    return p


if __name__ == "__main__" :
    """
    
    """
    args = parser().parse_args()

    if args.mode == "single":
    
        if args.mesh is None :
            raise ValueError("Use --mesh path/to/file.msh")
    
        solve_one(args.mesh, args.p, args.q, args.bx, args.by, args.c, args.out, args.plot)

    elif args.mode == "convergence" :

        if not args.meshes :
            raise ValueError("Use --meshes mesh1.msh mesh2.msh ...")
        
        run_convergence(args)
