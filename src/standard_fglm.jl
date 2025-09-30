@doc raw"""
    fglm_lex_shape_position(
        G::Vector{T},
        B::Vector{T}
    ) where T <: MPolyRingElem

Computes the lexicographic Gröbner basis of ideal ``(G) \subset K [x_1, \dots, x_n]`` via the standard FGLM algorithm.
If the lexicographic Gröbner basis of ``(G)`` is not in ``x_n``-shape position an error is returned.

# Arguments
- `G::Vector{MPolyRingElem}`: zero-dimensional Gröbner basis for which one wants to compute the lexicographic Gröbner basis.
- `B::Vector{MPolyRingElem}`: vector space basis of the quotient ring.

# Examples

```jldoctest
julia> P, (x, y) = polynomial_ring(GF(101), ["x", "y"], internal_ordering=:degrevlex);

julia> G = [y^2 + x - 1 + x, x^2 - y + 1];

julia> B = [P(1), x,  y, x * y];

julia> fglm_lex_shape_position(G, B)
Gröbner basis with elements
  1: y^4 + 99*y^2 + 97*y + 5
  2: x + 51*y^2 + 50
with respect to the ordering
  lex([x, y])
```
"""
function fglm_lex_shape_position(G::Vector{T}, B::Vector{T}) where T <: MPolyRingElem
    P = parent(first(G))
    R = base_ring(P)

    isa(R, Oscar.AbstractAlgebra.Field) || error("The FGLM algorithm requires a coefficient ring that is a field.")
    is_global(default_ordering(P)) || error("Default ordering of polynomial ring must be global.")
    is_zero_dimensional_gb(G) || error("Groebner basis is not zero-dimensional.")

    variables = gens(P)
    M_update = sparse_matrix(R)
    M_original = sparse_matrix(R)
    x_last = last(gens(P))
    f = one(P)
    d::Int = 0
    row_f = sparse_row(R, [(findfirst(==(one(P)), B), one(R))])
    push!(M_update, row_f)
    push!(M_original, row_f)
    r_old, M_update = rref(M_update)

    linearly_dependent = false
    while !linearly_dependent
        f *= x_last
        d += 1
        f = divrem(f, G)[2]
        row_f = normal_form_to_sparse_row(f, B)
        push!(M_update, row_f)

        r, M_update = rref(M_update)
        if r == r_old
            linearly_dependent = true
        else
            r_old = r
            push!(M_original, row_f)
        end
    end

    gb_lex = Vector{T}()
    x = solve(M_original, row_f)
    f = x_last^d - sum(el -> x_last^(el[1] - 1) * el[2], x)
    total_degree(f) == length(B) || error("LEX Groebner basis is not in shape position with respect to last variable.")
    push!(gb_lex, f)
    n = length(variables)
    for i in 1:(n - 1)
        row_f = normal_form_to_sparse_row(divrem(variables[i], G)[2], B)
        x = solve(M_original, row_f)
        f = variables[i] - sum(el -> x_last^(el[1] - 1) * el[2], x)
        push!(gb_lex, f)
    end

    return IdealGens(gb_lex, lex(P); isGB=true)
end

@doc raw"""
    fglm_lex_shape_position(
        I::MPolyIdeal
    )

Computes the lexicographic Gröbner basis of an ideal ``I \subset K [x_1, \dots, x_n]`` via the standard FGLM algorithm.
If the lexicographic Gröbner basis of ``I`` is not in ``x_n``-shape position an error is returned.

# Arguments
- `I::MPolyIdeal`: zero-dimensional ideal for which one wants to compute the lexicographic Gröbner basis.

# Examples

```jldoctest
julia> P, (x, y) = polynomial_ring(GF(101), ["x", "y"], internal_ordering=:degrevlex);

julia> I = ideal([y^2 + x - 1 + x, x^2 - y + 1]);

julia> fglm_lex_shape_position(I)
Gröbner basis with elements
  1: y^4 + 99*y^2 + 97*y + 5
  2: x + 51*y^2 + 50
with respect to the ordering
  lex([x, y])
```
"""
function fglm_lex_shape_position(I::MPolyIdeal)
    P = base_ring(I)
    R = base_ring(P)

    isa(R, Oscar.AbstractAlgebra.Field) || error("The FGLM algorithm requires a coefficient ring that is a field.")
    is_global(default_ordering(P)) || error("Default ordering of polynomial ring must be global.")

    if !haskey(I.gb, default_ordering(P)) || I.gb[default_ordering(P)].isReduced == false
        standard_basis(I, complete_reduction=true)
    end

    A, _ = quo(P, I)
    B = monomial_basis(A)

    I.gb[lex(P)] = fglm_lex_shape_position(gens(I.gb[default_ordering(P)]), B)

    return I.gb[lex(P)]
end
