using Oscar
using FglmAlgorithms
using Test
using Documenter

@testset verbose=true "All tests" begin
  doctestsetup = quote
    using Oscar
    using FglmAlgorithms
  end
  DocMeta.setdocmeta!(FglmAlgorithms, :DocTestSetup, doctestsetup; recursive=true)
  doctest(FglmAlgorithms; manual = false)
  set_verbosity_level(:groebner_walk, 0)

  include("fglm_standard.jl")
  include("fglm_sparse.jl")
  include("fglm_matrix.jl")
end

nothing
