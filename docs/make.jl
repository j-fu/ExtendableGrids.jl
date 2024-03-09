using Documenter, ExtendableGrids, ExampleJuggler, Gmsh
import CairoMakie
CairoMakie.activate!(; type = "svg", visible = false)
ExampleJuggler.verbose!(true)
ExtendableGridsGmshExt=Base.get_extension(ExtendableGrids, :ExtendableGridsGmshExt)
function mkdocs()
    cleanexamples()
    generated_examples = @docscripts(joinpath(@__DIR__, "..", "examples"),
                                     ["examples1d.jl", "examples2d.jl", "examples3d.jl", "gmsh.jl"],
                                     Plotter=CairoMakie)
    makedocs(; sitename = "ExtendableGrids.jl",
             modules = [ExtendableGrids,ExtendableGridsGmshExt],
             clean = false,
             doctest = true,
             authors = "J. Fuhrmann, Ch. Merdon, J. Taraz",
             repo = "https://github.com/j-fu/ExtendableGrids.jl",
             pages = [
                 "Home" => "index.md",
                 "extendablegrid.md",
                 "adjacency.md",
                 "vectorofconstants.md",
                 "typehierarchy.md",
                 "elementgeometry.md",
                 "shape_specs.md",
                 "coordinatesystem.md",
                 "subgrid.md",
                 "more.md",
                 "voronoi.md",
                 "assembly.md",
                 "cellfinder.md",
                 "arraytools.md",
                 "gridconstructors.md",
                 "output.md",
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
