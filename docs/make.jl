ENV["MPLBACKEND"] = "agg"
using Documenter, ExtendableGrids, Literate, GridVisualize, SimplexGridFactory, Gmsh
import CairoMakie, Triangulate

CairoMakie.activate!(; type = "svg", visible = false)

example_md_dir = joinpath(@__DIR__, "src", "examples")

examples1d = joinpath(@__DIR__, "..", "examples", "examples1d.jl")
include(examples1d)
examples2d = joinpath(@__DIR__, "..", "examples", "examples2d.jl")
include(examples2d)
examples3d = joinpath(@__DIR__, "..", "examples", "examples3d.jl")
include(examples3d)

include("makeplots.jl")

function mkdocs()
    Literate.markdown(examples1d, example_md_dir; documenter = false, info = false)
    Literate.markdown(examples2d, example_md_dir; documenter = false, info = false)
    Literate.markdown(examples3d, example_md_dir; documenter = false, info = false)

    generated_examples = joinpath.("examples", filter(x -> endswith(x, ".md"), readdir(example_md_dir)))

    makeplots(example_md_dir; Plotter = CairoMakie)

    makedocs(; sitename = "ExtendableGrids.jl",
             modules = [ExtendableGrids, Base.get_extension(ExtendableGrids, :ExtendableGridsGmshExt)],
             clean = false,
             warnonly = true,
             doctest = true,
             authors = "J. Fuhrmann, Ch. Merdon",
             repo = "https://github.com/j-fu/ExtendableGrids.jl",
             pages = [
                 "Home" => "index.md",
                 "adjacency.md",
                 "vectorofconstants.md",
                 "typehierarchy.md",
                 "elementgeometry.md",
                 "shape_specs.md",
                 "coordinatesystem.md",
                 "extendablegrid.md",
                 "subgrid.md",
                 "more.md",
                 "voronoi.md",
                 "assembly.md",
                 "cellfinder.md",
                 "arraytools.md",
                 "gridconstructors.md",
                 "refinement.md",
                 "regionedit.md",
                 "tokenstream.md",
                 "gmsh.md",
                 "allindex.md",
                 "Examples" => generated_examples,
             ])
end

mkdocs()

deploydocs(; repo = "github.com/j-fu/ExtendableGrids.jl.git")
