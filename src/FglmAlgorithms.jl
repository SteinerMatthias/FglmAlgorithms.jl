module FglmAlgorithms

using Oscar
import Oscar:
    base_ring,
    characteristic,
    coeff,
    coefficient_ring,
    divides,
    divrem,
    hcat,
    IdealGens,
    leading_monomial,
    matrix,
    monomial_basis
    MPolyRingElem,
    MPolyIdeal,
    monomials,
    one,
    parent,
    quo,
    rref,
    sparse_matrix,
    sparse_row,
    solve,
    standard_basis,
    total_degree,
    zero

include("berlekamp_massey.jl")
include("common.jl")
include("multiplication_matrix.jl")
include("standard_fglm.jl")
include("matrix_fglm.jl")
include("sparse_fglm.jl")

export fglm_lex_shape_position
export fglm_matrix_lex_shape_position
export fglm_sparse_lex_shape_position

end
