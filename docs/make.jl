push!(LOAD_PATH,"../src/")
using Documenter, ExtendableGrids

function mkdocs()
    makedocs(sitename="ExtendableGrids.jl",
             modules = [ExtendableGrids],
             doctest = true,
             clean = true,
             authors = "J. Fuhrmann, Ch. Merdon",
             repo="https://github.com/j-fu/ExtendableGrids.jl",
             pages=[
                 "Home"=>"index.md",
                 "adjacency.md",
                 "vectorofconstants.md",
                 "typehierarchy.md",
                 "elementgeometry.md",
                 "coordinatesystem.md",
                 "extendablegrid.md",
                 "subgrid.md",
                 "regionedit.md",
                 "simplexgrid.md",
                 "plot.md",
                 "tokenstream.md",
                 "allindex.md"
             ])
end

mkdocs()

deploydocs(repo = "github.com/j-fu/ExtendableGrids.jl.git")

