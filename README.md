# Numerical Experiments for ESDG and ESOFDG Methods

This repository contains the numerical experiment codes and mesh data accompanying the manuscript **"Error Estimates for Entropy-Stable Discontinuous Galerkin Methods for Hyperbolic Conservation Laws"** by Xu-Kun Chen, Yong Liu, and Chi-Wang Shu. The programs are intended to support the reproducibility of the reported experiments involving entropy-stable discontinuous Galerkin (ESDG) and entropy-stable oscillation-free discontinuous Galerkin (ESOFDG) methods.

## Repository Contents

| File | Description |
| --- | --- |
| `00_1d-scalar.cpp` | C++ implementation for one-dimensional scalar conservation-law experiments. |
| `00_2d-Euler-PBC.cpp` | C++ implementation for the two-dimensional Euler experiment with periodic boundary conditions (PBC). |
| `00_2d-Euler-version-1.cpp` | C++ implementation for the two-dimensional Euler experiment with exact-solution boundary data. |
| `00_CaculateDeter.py` | Auxiliary Python script for determinant calculation. |
| `square.zip` | Complete collection of the two-dimensional triangular mesh files used in the numerical experiments. |

## Requirements

- A C++ compiler supporting **C++17** (e.g., `g++`).
- Python 3 for the optional determinant utility.
- A utility for extracting `.zip` archives.

## Setup and Compilation

Download or clone the repository, then extract `square.zip` so that the `.msh` files are located in the working directory from which the two-dimensional programs will be run.

Compile the desired C++ source file using a C++17 compiler. For example:

```bash
g++ -std=c++17 -O2 00_2d-Euler-PBC.cpp -o euler_pbc
```

Then run the corresponding executable.

To compute the determinants of the Vandermonde matrices for $k=3,4$, run the auxiliary Python script:

```bash
python3 00_CaculateDeter.py
```

These calculations provide a computational verification of the nonsingularity of the corresponding Vandermonde matrices, supporting the proposition in the manuscript.

## Reference

Xu-Kun Chen, Yong Liu, and Chi-Wang Shu, *Error Estimates for Entropy-Stable Discontinuous Galerkin Methods for Hyperbolic Conservation Laws* (manuscript).
