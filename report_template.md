
# Report template — FEM Project

## 1. Problem

We study

\[
-\Delta u + b\cdot \nabla u + cu = f
\quad \text{in } \Omega,
\qquad
u=0 \quad \text{on } \partial\Omega.
\]

The approximation space is the conforming P1 finite element space on a triangular mesh.

---

## 2. Variational formulation

Let \(V=H_0^1(\Omega)\). Multiply the PDE by \(v\in V\) and integrate:

\[
\int_\Omega (-\Delta u)v
+
\int_\Omega (b\cdot\nabla u)v
+
c\int_\Omega uv
=
\int_\Omega fv.
\]

By integration by parts,

\[
\int_\Omega (-\Delta u)v
=
\int_\Omega \nabla u\cdot\nabla v
-
\int_{\partial\Omega} \partial_nu\,v.
\]

Since \(v=0\) on \(\partial\Omega\), the boundary term vanishes. Therefore:

\[
a(u,v)=\ell(v)
\]

where

\[
a(u,v)
=
\int_\Omega \nabla u\cdot\nabla v
+
\int_\Omega (b\cdot\nabla u)v
+
c\int_\Omega uv,
\]

and

\[
\ell(v)=\int_\Omega fv.
\]

---

## 3. Discretization

Let

\[
u_h=\sum_j U_j\phi_j.
\]

The linear system is

\[
(K+C+cM)U=F.
\]

The matrices are

\[
M_{ij}=\int_\Omega \phi_j\phi_i,
\]

\[
K_{ij}=\int_\Omega \nabla\phi_j\cdot\nabla\phi_i,
\]

\[
C_{ij}=\int_\Omega (b\cdot\nabla\phi_j)\phi_i.
\]

---

## 4. Elementary convection matrix

On one triangle \(T\),

\[
(C_T)_{ij}
=
\int_T (b\cdot\nabla\phi_j)\phi_i.
\]

For P1 elements, \(\nabla\phi_j\) is constant on \(T\), so

\[
(C_T)_{ij}
=
(b\cdot\nabla\phi_j)\int_T \phi_i.
\]

Since

\[
\int_T \phi_i = \frac{|T|}{3},
\]

we obtain

\[
(C_T)_{ij}
=
\frac{|T|}{3}(b\cdot\nabla\phi_j).
\]

This is the formula implemented in `elementary_convection_matrix`.

---

## 5. Exact solution

The source is

\[
f(x,y)=e^{(b_xx+b_yy)/2}\sin(p\pi x)\sin(q\pi y).
\]

We look for

\[
u_{ex}=\alpha f.
\]

One obtains

\[
(-\Delta + b\cdot\nabla + c)f
=
\left[
\pi^2(p^2+q^2)+c+\frac{|b|^2}{4}
\right]f.
\]

Therefore,

\[
\alpha =
\frac{1}{
\pi^2(p^2+q^2)+c+\frac{|b|^2}{4}
}.
\]

---

## 6. Error formulas

The project asks for

\[
\frac{\|u_h-\Pi_hu_{ex}\|_{L^2}}
{\|\Pi_hu_{ex}\|_{L^2}}
\]

and

\[
\frac{\|u_h-\Pi_hu_{ex}\|_{H^1}}
{\|\Pi_hu_{ex}\|_{H^1}}.
\]

With vectors \(U_h\) and \(U_{ex}\),

\[
\|u_h-\Pi_hu_{ex}\|_{L^2}^2
=
(U_h-U_{ex})^T M (U_h-U_{ex}),
\]

and

\[
\|u_h-\Pi_hu_{ex}\|_{H^1}^2
=
(U_h-U_{ex})^T(M+K)(U_h-U_{ex}).
\]

Expected rates:

\[
L^2 = O(h^2),
\qquad
H^1 = O(h).
\]

---

## 7. Figures to include

- Mesh plot.
- Numerical solution \(u_h\).
- Exact solution sampled at nodes.
- Log-log convergence plot.
- Boundary connected components for `barwith4holes.msh`.

---

## 8. Connected components

Boundary edges form a graph. Two edges belong to the same component if one can move from one to the other through shared boundary vertices. A BFS/DFS algorithm labels each component.

This is used to impose

\[
u=j \quad \text{on } \Gamma_j.
\]

---

## 9. Discussion

Increasing \(|b|\) makes the problem convection-dominated. Standard Galerkin P1 FEM can then show oscillations on coarse meshes. This is a classical motivation for stabilized methods such as SUPG.
