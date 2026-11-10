import numpy as np
from itertools import permutations

def triangle_nodes(k):
    z, _ = np.polynomial.legendre.leggauss(k+1)
    z = (z+1)/2

    nodes = []
    for t in z:
        nodes.extend([
            (0, t, 1-t),
            (t, 0, 1-t),
            (t, 1-t, 0)
        ])

    if k == 3:
        a = 0.1870738791912771
        b = 0.5841571139756569
        abc = (a, b, 1-a-b)

        for p in permutations(abc):
            nodes.append(p)

    elif k == 4:
        for s in (0.4384239524408185,
                  0.1394337314154536):
            u = 1-2*s
            nodes.extend([
                (s, s, u),
                (s, u, s),
                (u, s, s)
            ])

        nodes.append((1/3, 1/3, 1/3))

    return np.array(nodes)

def polynomial_powers(k):
    return [
        (a, degree-a)
        for degree in range(k+1)
        for a in range(degree, -1, -1)
    ]

def vandermonde(k):
    if k == 3:
        powers = polynomial_powers(3) + [
            (0,4), (4,1), (1,4), (6,0),
            (5,1), (4,2), (2,4), (0,6)
        ]
    else:
        powers = polynomial_powers(5) + [(6,0)]

    nodes = triangle_nodes(k)

    X = 3*nodes[:, 2]-1
    Y = 3*nodes[:, 1]-1

    V = np.column_stack([
        X**a * Y**b for a,b in powers
    ])

    return V

for k in (3,4):
    V = vandermonde(k)

    print(f"k = {k}")
    print(f"Matrix size = {V.shape}")
    print(f"det(V) = {np.linalg.det(V):.12e}")
    print(f"cond(V) = {np.linalg.cond(V):.6e}")
    print()