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
Pages = ["multiplication-matrix.md"]
```

# Multiplication Matrix
Let $K$ be a field, let $I \subset K [x_1, \dots, x_n]$ be a zero-dimensional ideal, and let $G \subset I$ be a $>$-Gröbner basis of $I$.
Then $A = K [x_1, \dots, x_n] / I$ is a finite dimensional $K$-vector space, and a vector space basis is given by all monomials $B$ not contained in the ideal of $>$-leading terms of $I$.
In particular, for $f \in K [x_1, \dots, x_n]$ the map $x \mapsto x \cdot f$ is $K$-linear on $A$.
So, it can be represented by a matrix $\mathbf{M}_f \in K^{|B| \times |B|}$.
The following functions provide functionality for the computation of $\mathbf{M}_f$, given as input
- a polynomial $f$,
- a zero-dimensional $>$-Gröbner basis $G$, and
- the $K$-vector space basis $B$ of $K [x_1, \dots, x_n] / (G)$ induced by $G$.
For maximum performance, no boundary checks are performed.

## Functionality
```@docs
    normal_form_to_sparse_row(
        f::T,
        B::Vector{T}
    ) where T <: MPolyRingElem

    multiplication_matrix_dense(
        f::T,
        G::Vector{T},
        B::Vector{T}
    ) where T <: MPolyRingElem

    multiplication_matrix(
        f::T,
        G::Vector{T},
        B::Vector{T}
    ) where T <: MPolyRingElem
```
