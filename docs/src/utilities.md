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
Pages = ["utilities.md"]
```
# Utilities
Collection of functions which are shared by the FGLM algorithms.

## Functionality
```@docs
    is_zero_dimensional_gb(
        G::Vector{T}
    ) where T <: MPolyRingElem
```
