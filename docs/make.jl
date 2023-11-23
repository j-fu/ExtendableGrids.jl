using Documenter, ExtendableGrids, ExampleJuggler, Gmsh
import CairoMakie
CairoMakie.activate!(; type = "svg", visible = false)
ExampleJuggler.verbose!(true)

function mkdocs()
    cleanexamples()
    generated_examples = @docscripts(joinpath(@__DIR__, "..", "examples"),
                                     ["examples1d.jl", "examples2d.jl", "examples3d.jl", "gmsh.jl"],
                                     Plotter=CairoMakie)
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
