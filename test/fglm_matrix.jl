@testset "Matrix FGLM" begin
    K = GF(101)
    P, (x, y, z) = polynomial_ring(K, ["x", "y", "z"], internal_ordering=:degrevlex)

    I = ideal([x^2 - x * y + 5 * z^2 - x + 2 * y - 3, x * y - 3, z^2 + x + z + y + 5])
    gb_lex = fglm_matrix_lex_shape_position(I)

    gb_lex_expected = [z^6 + 35*z^5 + 44*z^4 + 90*z^3 + 24*z^2 + 67*z + 32,
                       y + 18*z^5 + 51*z^4 + 51*z^3 + 11*z^2 + 25*z + 69,
                       x + 83*z^5 + 50*z^4 + 50*z^3 + 91*z^2 + 77*z + 37]
    isGB = all(gens(gb_lex)) do f
      f in gb_lex_expected
    end

    @test isGB==true
end
