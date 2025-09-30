@doc raw"""
    berlekamp_massey(
        a::Vector{S},
        x::T
    ) where {S <: RingElement, T <: MPolyRingElem}

Computes the minimal polynomial of a linear recurrent sequence.
**For internal use only.**

# Arguments
- `a::Vector{S}`: a sequence of ring (in practice only field) elements.
- `x::T`: variable in which the polynomial is computed.
"""
function berlekamp_massey(a::Vector{S}, x::T) where {S <: RingElement, T <: MPolyRingElem}
    R = base_ring(parent(x))
    m = length(a) << -1

    f_0 = zero(R)
    for el in reverse(a)
        f_0 *= x
        f_0 += el
    end
    f_1 = x^(2 * m)
    s_0 = one(R)
    s_1 = zero(R)

    while total_degree(f_1) >= m
        q, r = divrem(f_0, f_1)
        f_0 = f_1
        f_1 = r
        s_0, s_1 = (s_1, s_0 - q * s_1)
    end

    f = zero(R)
    for (mon, coeff) in zip(monomials(s_1), reverse(collect(coefficients(s_1))))
        f += mon * coeff
    end

    return f / leading_coefficient(f)
end
