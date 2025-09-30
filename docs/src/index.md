```@meta
CurrentModule = FglmAlgorithms
```

# Getting Started

[FglmAlgorithms.jl](https://github.com/SteinerMatthias/FglmAlgorithms.jl) implements term order conversion algorithms for zero-dimensional ideals to a lexicographic Gröbner basis on top of the [OSCAR](https://www.oscar-system.org/) computer algebra system.

## Installation
The package requires [Julia](https://julialang.org/) 1.11 or higher and [OSCAR](https://www.oscar-system.org/) 1.0 or higher.

For installation, simply type
```julia
julia> using Pkg; Pkg.add("FglmAlgorithms")
```

## A Simple Example
After installation you can use the package in combination with [OSCAR](https://www.oscar-system.org/).
We provide a simple example
```julia
using Oscar, FglmAlgorithms

P, (x, y) = polynomial_ring(GF(101), ["x", "y"], internal_ordering=:degrevlex)
I = ideal([y^2 - x + 1, x^2 + x - y - 1])

gb_lex_1 = fglm_lex_shape_position(I)
gb_lex_2 = fglm_matrix_lex_shape_position(I)
gb_lex_3 = fglm_sparse_lex_shape_position(I)
```