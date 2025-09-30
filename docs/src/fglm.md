```@meta
CurrentModule = FglmAlgorithms
DocTestSetup = quote
  using FglmAlgorithms
end
```

```@setup fglm
using FglmAlgorithms
```

```@contents
Pages = ["fglm.md"]
```

# FGLM Algorithms
Let $K$ be a field, and let $I \subset K [x_1, \dots, x_n]$ be a zero-dimensional ideal.
Givean a $>_1$-Gröbner basis $G \subset I$ be a  of $I$, with variants of the FGLM algorithm one can perform term order conversion from $G$ to a $>_2$-Gröbner basis.
This module implements the following variants of FGLM:
- The original FGLM algorithm [https://doi.org/10.1006/jsco.1993.1051](https://doi.org/10.1006/jsco.1993.1051).
- FGLM with fast matrix multiplication [https://doi.org/10.1145/2608628.2608669](https://doi.org/10.1145/2608628.2608669).
- Sparse FGLM [https://doi.org/10.1016/j.jsc.2016.07.025](https://doi.org/10.1016/j.jsc.2016.07.025).

Currently, only the following FGLM algorithm functionality is implemented:
- Only term order conversion to a lexicographic order is supported.
- Only the probabilistic versions of the fast matrix multiplication and sparse FGLM are implemented.
- Fast matrix multiplication and sparse FGLM are only supported over finite fields.
- FGLM algorithms are only supported for ideals in $x_n$-shape position, i.e. the lexicographic Gröbner basis has the shape
```math
\begin{aligned}
    g_1 &= x_1 - f_1 (x_n), \\
    &\vdots \\
    g_{n - 1} &= x_{n - 1} f_{n - 1} (x_n), \\
    g_n &= f_n (x_n).
\end{aligned}
```

## Functionality

### Standard FGLM
```@docs
    fglm_lex_shape_position(
        I::MPolyIdeal
    )

    fglm_lex_shape_position(
        G::Vector{T},
        B::Vector{T}
    ) where T <: MPolyRingElem
```

### FGLM With Fast Matrix Multiplication
```@docs
    fglm_matrix_lex_shape_position(
        I::MPolyIdeal
    )

    fglm_matrix_lex_shape_position(
        G::Vector{T},
        B::Vector{T}
    ) where T <: MPolyRingElem
```

### Sparse FGLM
```@docs
    fglm_sparse_lex_shape_position(
        I::MPolyIdeal
    )

    fglm_sparse_lex_shape_position(
        G::Vector{T},
        B::Vector{T}
    ) where T <: MPolyRingElem
```
