using Oscar
using FglmAlgorithms
using Documenter

DocMeta.setdocmeta!(FglmAlgorithms,
                    :DocTestSetup,
                    :(using Oscar, FglmAlgorithms);
                    recursive=true,
                   )

makedocs(;
  modules=[FglmAlgorithms],
  doctest=true,
  linkcheck=true,
  checkdocs=:exports,
  authors="SteinerMatthias",
  sitename="FglmAlgorithms.jl",
  format=Documenter.HTML(;prettyurls=true),
  pages = [
    "Home" => "index.md",
    "Functionality" => ["fglm.md",
                        "multiplication-matrix.md",
                        "berlekamp-massey.md",
                        "utilities.md",
                       ],
  ],
)

deploydocs(
  repo="https://github.com/SteinerMatthias/FglmAlgorithms.jl",
)
