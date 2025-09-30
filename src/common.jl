@doc raw"""
    is_zero_dimensional_gb(
        G::Vector{T}
    ) where T <: MPolyRingElem

A fast check whether a polynomial vector is a zero-dimensional Gröbner basis with respect to the default ordering of the polynomial ring.
**For internal use only.**

# Arguments
- `G::Vector{T}`: an asserted zero-dimensional Gröbner basis.
"""
function is_zero_dimensional_gb(G::Vector{T}) where T <: MPolyRingElem
    P = parent(first(G))
    variables = gens(P)
    lms = map(leading_monomial, G)
    lms = filter(mon -> is_univariate(mon), lms)
    for var in variables
        var_divides = false
        for mon in lms
            var_divides = divides(mon, var)[1]
            if var_divides
                break
            end
        end
        if !var_divides
            return false
        end
    end
    return true
end
