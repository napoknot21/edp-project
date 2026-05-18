# FEM Project

Sorbonne Université project for solving a 2D finite-element test problem on triangular
meshes. The code assembles P1 mass, stiffness and convection matrices, applies
Dirichlet boundary conditions, compares the numerical solution with an exact
solution and can generate convergence plots.

## Project Layout

```
project/
  main.py                  # command-line entry point
  requirements.txt         # Python dependencies
  src/
    mesh_io.py             # mesh loading and boundary extraction
    geometry.py            # triangle/edge geometry helpers
    matrices.py            # local and global FEM matrices
    load.py                # load-vector assembly
    dirichlet.py           # Dirichlet boundary solvers
    exact.py               # source and exact solution
    errors.py              # mesh size and relative error norms
    plotting.py            # PNG plots
  data/                    # example .msh meshes and .geo file
  outputs/                 # generated figures
```

## Installation

From the `project` directory:

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

On macOS/Linux, activate the environment with:

```bash
source .venv/bin/activate
```

## Run One Mesh

Basic command:

```bash
python main.py --mesh data/Meshes-20260518/Lshape.msh
```

With all main parameters and PNG plots:

```bash
python main.py ^
  --mode single ^
  --mesh data/Meshes-20260518/disk1hole.msh ^
  --p 2 ^
  --q 3 ^
  --bx 1 ^
  --by 1 ^
  --c 5 ^
  --out outputs ^
  --plot
```

The command prints the mesh size `h`, the coefficient `alpha`, and the relative
L2/H1 errors. With `--plot`, it also saves:

- `mesh_<mesh_name>.png`
- `uh_<mesh_name>.png`
- `uex_nodes_<mesh_name>.png`

## Convergence Mode

Use several meshes to compare the error as the mesh size decreases:

```bash
python main.py ^
  --mode convergence ^
  --meshes data/square_h0.2.msh data/square_h0.1.msh data/square_h0.05.msh data/square_h0.025.msh ^
  --p 1 ^
  --q 1 ^
  --bx 1 ^
  --by 1 ^
  --c 1 ^
```

This saves a figure named like:

```text
outputs/convergence_p1_q1_bx1_by1_c1.png
```

## Command-Line Options

```text
--mode      single or convergence, default: single
--mesh      path to one .msh file, used in single mode
--meshes    list of .msh files, used in convergence mode
--p         integer index in the exact/source function, default: 1
--q         integer index in the exact/source function, default: 1
--bx        x-component of the convection vector, default: 1.0
--by        y-component of the convection vector, default: 1.0
--c         reaction coefficient, default: 1.0
--out       output directory for generated PNG files, default: outputs
--plot      save mesh and solution figures in single mode
```

You can also ask Python for the generated CLI help:

```bash
python main.py --help
```

## Pydoc

The code contains short docstrings for the main public functions. To inspect
them in the terminal:

```bash
python -m pydoc main
python -m pydoc src.matrices
python -m pydoc src.dirichlet
```

To generate HTML documentation files:

```bash
python -m pydoc -w main src.mesh_io src.geometry src.matrices src.load src.dirichlet src.exact src.errors src.plotting
```

This creates `.html` files in the current directory. They can be opened in a
browser.

## Mesh Examples

Included meshes:

```text
data/square_h0.2.msh
data/square_h0.1.msh
data/square_h0.05.msh
data/square_h0.025.msh
data/Meshes-20260518/Lshape.msh
data/Meshes-20260518/disk1hole.msh
data/Meshes-20260518/fish.msh
data/Meshes-20260518/barwith4holes.msh
data/Meshes-20260518/rectangle01.msh
```

