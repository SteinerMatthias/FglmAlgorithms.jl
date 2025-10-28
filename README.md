# FglmAlgorithms.jl
FglmAlgorithms.jl is a Julia package providing implementations of FGLM [[1](#1), [2](#2), [3](#3)] term order conversion algorithms on top of [OSCAR](https://www.oscar-system.org/).
Currently, only term order conversion to **lexicographic** Gröbner bases in **shape position** is supported.
I.e., the lexicographic Gröbner basis must be of the form
```math
\begin{aligned}
    f_1 &= x_1 - g_1 (x_n), \\
    &\vdots \\
    f_{n - 1} &= x_{n - 1} - g_{n - 1} (x_n), \\
    f_n &= g_n (x_n).
\end{aligned}
```
Moreover, for the algorithms of [[2](#2), [3](#3)] **only probabilistic variants** are implemented.

## Installation
FglmAlgorithms.jl can be installed from the Julia REPL via
```julia
using Pkg; Pkg.add(url="https://github.com/SteinerMatthias/FglmAlgorithms.jl")
```

## Usage
FglmAlgorithmsGLM.jl exports three functions: `fglm_lex_shape_position`, `fglm_matrix_lex_shape_position` and `fglm_sparse_lex_shape_position`.
- `fglm_lex_shape_position`: Term order conversion via the standard FGLM algorithm [[1](#1)].
- `fglm_matrix_lex_shape_position`: Term order conversion via fast matrix multiplication [[2](#2)]. This function handles multiplicaton matrices as dense matrices.
- `fglm_sparse_lex_shape_position`: Term order conversion via the sparse FGLM algorithm [[3](#3)]. This function handles multiplicaton matrices via [OSCAR](https://www.oscar-system.org/)'s sparse matrix interface.

As input the FGLM functions either require
- a zero-dimensional ideal $I \subset P$ (the $>$-Gr\"obner basis of $I$ is first computed if it is not known), or
- the $>$-Gr\"obner basis $G \subset I$ and the $K$-vector space basis $B$ of $P / (G)$ as vectors of polynomials.

First, one needs to initalize an ideal.
```julia
using Oscar

K = GF(101)
P, (x, y, z) = polynomial_ring(K, ["x", "y", "z"], internal_ordering=:degrevlex)
G = [y^2 + y * z - z^2 + x, x^2 + 50 * z^2 + 1, z^2 + x - y] # This is already a DRL Gröbner basis
I = ideal(G)
# A vectos space basis of I is given by
B = [P(1), x, y, z, x * y, x * z, y * z, x * y * z]
```
Now we compute the lexicographic Gröbner basis.
```julia
using FglmAlgorithms

G_lex = fglm_lex_shape_position(I)
G_lex = fglm_matrix_lex_shape_position(I)
G_lex = fglm_sparse_lex_shape_position(I)
```
Alternatively, we can compute it with $G$ and $B$ as inputs.
```julia
G_lex = fglm_lex_shape_position(G, B)
G_lex = fglm_matrix_lex_shape_position(G, B)
G_lex = fglm_sparse_lex_shape_position(G, B)
```

## Documentation

The documentation for the FglmAlgorithms.jl package can be found at [https://steinermatthias.github.io/FglmAlgorithms.jl/](https://steinermatthias.github.io/FglmAlgorithms.jl/).

## References
<a id="1">[1]</a>
J.C. Faugère, P. Gianni, D. Lazard, T. Mora.
(1993).
Efficient Computation of Zero-dimensional Gröbner Bases by Change of Ordering.
Journal of Symbolic Computation.
Volume 16, Issue 4.
https://doi.org/10.1006/jsco.1993.1051.

<a id="2">[2]</a>
J.C. Faugère, P. Gaudry, L, Huot, G. Renault.
(2014).
Sub-cubic change of ordering for Gröbner basis: a probabilistic approach.
Proceedings of the 39th International Symposium on Symbolic and Algebraic Computation (ISSAC '14)
https://doi.org/10.1145/2608628.2608669.

<a id="3">[3]</a>
J.C. Faugère, C. Mou.
(2017).
Sparse FGLM algorithms.
Journal of Symbolic Computation.
Volume 80, Part 3.
https://doi.org/10.1016/j.jsc.2016.07.025.
