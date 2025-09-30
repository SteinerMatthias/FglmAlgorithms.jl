@doc raw"""
    fglm_lex_shape_fglm_sparse_lex_shape_positionposition(
        G::Vector{T},
        B::Vector{T}
    ) where T <: MPolyRingElem

Computes the lexicographic Gröbner basis of ideal ``(G) \subset K [x_1, \dots, x_n]`` via the sparse FGLM algorithm.
If the lexicographic Gröbner basis of ``(G)`` is not in ``x_n``-shape position an error is returned.

# Arguments
- `G::Vector{MPolyRingElem}`: zero-dimensional Gröbner basis for which one wants to compute the lexicographic Gröbner basis.
- `B::Vector{MPolyRingElem}`: vector space basis of the quotient ring.

# Examples

```jldoctest
julia> P, (x, y) = polynomial_ring(GF(101), ["x", "y"], internal_ordering=:degrevlex);

julia> G = [y^2 + x - 1 + x, x^2 - y + 1];

julia> B = [P(1), x,  y, x * y];

julia> fglm_sparse_lex_shape_position(G, B)
Gröbner basis with elements
  1: x + 51*y^2 + 50
  2: y^4 + 99*y^3 + 97*y^2 + 5
with respect to the ordering
  lex([x, y])
```
"""
function fglm_sparse_lex_shape_position(G::Vector{T}, B::Vector{T}) where T <: MPolyRingElem
    P = parent(first(G))
    R = base_ring(P)

    isa(R, Oscar.AbstractAlgebra.Field) || error("The FGLM algorithm requires a coefficient ring that is a field.")
    characteristic(R) > 0 || error("Only finite fields are supported.")
    is_global(default_ordering(P)) || error("Default ordering of polynomial ring must be global.")
    is_zero_dimensional_gb(G) || error("Groebner basis is not zero-dimensional.")

    variables = gens(P)
    n = length(variables)
    R = base_ring(P)
    x_last = last(gens(P))
    D = length(B)
    index_one = findfirst(==(one(P)), B)
    indices_variables = map(i -> findfirst(==(variables[i]), B), 1:(n - 1))

    M = multiplication_matrix(x_last, G, B)
    r = [rand(R) for _ in 1:D]
    s = [r[index_one]]


    B_i = [Vector{typeof(zero(R))}() for _ in 1:n - 1]
    for i in 1:(n - 1)
        push!(B_i[i], r[indices_variables[i]])
    end

    for i in 1:(2 * D - 1)
        r = M * r
        push!(s, r[index_one])
        if i < D
            for j in 1:(n - 1)
                push!(B_i[j], r[indices_variables[j]])
            end
        end
    end

    gb_lex = Vector{T}()

    H = [s[1:D]]
    for i in 2:D
        push!(H, s[i:D + i - 1])
    end
    H = matrix(H)
    #H = inv(matrix(H))

    for i in 1:(n - 1)
        b_i = B_i[i]
        c_i = solve(H, b_i)
        #c_i = H * b_i
        f_i = zero(P)
        for c in reverse(c_i)
            f_i *= x_last
            f_i += c
        end
        push!(gb_lex, variables[i] - f_i)
    end

    f = berlekamp_massey(s, x_last)
    total_degree(f) == D || error("LEX Groebner basis is not in shape position with respect to last variable.")
    push!(gb_lex, f)

    return Oscar.IdealGens(gb_lex, lex(P); isGB=true)
end

@doc raw"""
    fglm_sparse_lex_shape_position(
        I::MPolyIdeal
    )

Computes the lexicographic Gröbner basis of an ideal ``I \subset K [x_1, \dots, x_n]`` via the sparse FGLM algorithm.
If the lexicographic Gröbner basis of ``I`` is not in ``x_n``-shape position an error is returned.

# Arguments
- `I::MPolyIdeal`: zero-dimensional ideal for which one wants to compute the lexicographic Gröbner basis.

# Examples

```jldoctest
julia> P, (x, y) = polynomial_ring(GF(101), ["x", "y"], internal_ordering=:degrevlex);

julia> I = ideal([y^2 + x - 1 + x, x^2 - y + 1]);

julia> fglm_sparse_lex_shape_position(I)
Gröbner basis with elements
  1: x + 51*y^2 + 50
  2: y^4 + 99*y^3 + 97*y^2 + 5
with respect to the ordering
  lex([x, y])
```
"""
function fglm_sparse_lex_shape_position(I::MPolyIdeal)
    P = base_ring(I)
    R = base_ring(P)

    isa(R, Oscar.AbstractAlgebra.Field) || error("The FGLM algorithm requires a coefficient ring that is a field.")
    characteristic(R) > 0 || error("Only finite fields are supported.")
    is_global(default_ordering(P)) || error("Default ordering of polynomial ring must be global.")

    if !haskey(I.gb, default_ordering(P)) || I.gb[default_ordering(P)].isReduced == false
        standard_basis(I, complete_reduction=true)
    end

    A, _ = quo(P, I)
    B = monomial_basis(A)

    I.gb[lex(P)] = fglm_sparse_lex_shape_position(gens(I.gb[default_ordering(P)]), B)

    return I.gb[lex(P)]
end
