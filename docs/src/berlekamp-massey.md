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
Pages = ["berlekamp-massey.md"]
```

# Berlekamp-Massey Algorithm
The algorithm is used to compute the univariate polynomial in the fast matrix multiplication and sparse FGLM algorithms.

## Functionality
```@docs
    berlekamp_massey(
        a::Vector{S},
        x::T
    ) where {S <: RingElement, T <: MPolyRingElem}
```
