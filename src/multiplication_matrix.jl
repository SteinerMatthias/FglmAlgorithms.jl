@doc raw"""
    normal_form_to_sparse_row(
        f::T,
        B::Vector{T}
    ) where T <: MPolyRingElem

Computes the sparse coefficient vector of an element ``f`` in a zero-dimensional quotient ring ``K[x_1, ..., x_n] / I`` with respect the K-vector space basis ``B``.
Assumes that ``f`` is already in normal form.
**For internal use only.**

# Arguments
- `f::MPolyRingElem`: a multivariate polynomial in normal form.
- `B::Vector{MPolyRingElem}`: a vector space basis of ``K[x_1, ..., x_n] / I``.
"""
function normal_form_to_sparse_row(f::T, B::Vector{T}) where T <: MPolyRingElem
    R = base_ring(parent(f))
    row::Vector{Tuple{Int, RingElement}} = []
    for (i, b) in enumerate(B)
        if coeff(f, b) != zero(R)
            push!(row, (i, coeff(f, b)))
        end
    end
    return sparse_row(R, row)
end

@doc raw"""
    multiplication_matrix_sparse(
        f::T,
        G::Vector{T},
        B::Vector{T}
    ) where T <: MPolyRingElem

For a polynomial ``f``, computes the multiplication matrix of ``x \mapsto x \cdot f`` over ``K [x_1, ..., x_n] / (G)`` with respect to the K-vector space basis ``B``.
The multiplication matrix is represented as sparse matrix.
**For internal use only.**

# Arguments
- `f::MPolyRingElem`: a multivariate polynomial.
- `G::Vector{MPolyRingElem}`: a zero-dimensional Gröbner basis with respect to the default ordering of the polynomial ring.
- `B::Vector{MPolyRingElem}`: vector space basis of ``K[x_1, ..., x_n] / (G)`` with respect to the Gröbner basis G.
"""
function multiplication_matrix_sparse(f::T, G::Vector{T}, B::Vector{T}) where T <: MPolyRingElem
    R = base_ring(parent(f))
    M = sparse_matrix(R)
    for b in B
        r_fb = divrem(f * b, G)[2]
        row_fb = normal_form_to_sparse_row(r_fb, B)
        push!(M, row_fb)
    end
    return M
end

@doc raw"""
    multiplication_matrix_dense(
        f::T,
        G::Vector{T},
        B::Vector{T}
    ) where T <: MPolyRingElem

For a polynomial f, computes the multiplication matrix of ``x \mapsto x \cdot f`` over ``K [x_1, ..., x_n] / (G)`` with respect to the ``K``-vector space basis ``B``.
The multiplication matrix is represented as dense matrix.
**For internal use only**

# Arguments
- `f::MPolyRingElem`: a multivariate polynomial.
- `G::Vector{MPolyRingElem}`: a zero-dimensional Gröbner basis with respect to the default ordering of the polynomial ring.
- `B::Vector{MPolyRingElem}`: vector space basis of ``K [x_1, ..., x_n] / (G)`` with respect to the Gröbner basis ``G``.
"""
function multiplication_matrix_dense(f::T, G::Vector{T}, B::Vector{T}) where T <: MPolyRingElem
    R = base_ring(parent(f))
    M = Vector{Vector{RingElement}}()
    for b in B
        r_fb = divrem(f * b, G)[2]
        row_fb = map(c -> coeff(r_fb, c), B)
        push!(M, row_fb)
    end
    return matrix(M)
end

@doc raw"""
    multiplication_matrix(
        f::T,
        G::Vector{T},
        B::Vector{T};
        sparse::Bool=true
    ) where T <: MPolyRingElem

For a polynomial ``f``, computes the multiplication matrix of ``x \mapsto x \cdot f`` over ``K [x_1, ..., x_n] / (G)`` with respect to the ``K``-vector space basis ``B``.
**For internal use only.**

# Arguments
- `f::MPolyRingElem`: a multivariate polynomial.
- `G::Vector{MPolyRingElem}`: a zero-dimensional Gröbner basis with respect to the default ordering of the polynomial ring.
- `B::Vector{MPolyRingElem}`: vector space basis of ``K [x_1, ..., x_n] / (G)`` with respect to the Gröbner basis ``G``.
- `sparse::Bool`: optional parameter to reprsent the multiplication matrix as sparse matrix. Default valuse is true.
"""
function multiplication_matrix(f::T,
                               G::Vector{T},
                               B::Vector{T};
                               sparse::Bool=true) where T <: MPolyRingElem
    if sparse
        return multiplication_matrix_sparse(f, G, B)
    else
        return multiplication_matrix_dense(f, G, B)
    end
end
