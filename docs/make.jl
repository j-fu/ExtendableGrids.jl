push!(LOAD_PATH,"../src/")
using Documenter, XGrid

makedocs(sitename="XGrid.jl",
         modules = [XGrid],
         doctest = true,
         clean = true,
         authors = "J. Fuhrmann, Ch. Merdon",
         repo="https://github.com/j-fu/XGrid.jl",
         pages=[
             "Home"=>"index.md",
             "adjacency.md",
             "vectorofconstants.md",
             "elementgeometry.md",
             "coordinatesystem.md",
             "extendablegrid.md",
             "subgrid.md",
             "regionedit.md",
             "simplexgrid.md",
             "tokenstream.md",
             "plot.md",
             "allindex.md"
         ])

deploydocs(repo = "github.com/j-fu/XGrid.jl.git")

