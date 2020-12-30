ENV["MPLBACKEND"]="agg"
using Documenter, ExtendableGrids, Literate
import PyPlot


example_md_dir  = joinpath(@__DIR__,"src","examples")

examples1d=joinpath(@__DIR__,"..","examples","examples1d.jl")
include(examples1d)
examples2d=joinpath(@__DIR__,"..","examples","examples2d.jl")
include(examples2d)
examples3d=joinpath(@__DIR__,"..","examples","examples3d.jl")
include(examples3d)
plotting=joinpath(@__DIR__,"..","examples","plotting.jl")
include(plotting)




function makeplots(picdir)
    function makeplot(func)
        PyPlot.clf()
        f=getfield(Main,Symbol(func))
        plot(f(), Plotter=PyPlot)
        PyPlot.savefig(joinpath(picdir,func*".svg"))
    end

    function makeplot2(func)
        PyPlot.clf()
        f=getfield(Main,Symbol(func))
        f(Plotter=PyPlot)
        PyPlot.savefig(joinpath(picdir,func*".svg"))
    end

    makeplot("interval_from_vector")
    makeplot("interval_localref")
    makeplot("interval_multiregion")
    makeplot("interval_subgrid")
    makeplot("rectangle")
    makeplot("rectangle_localref")
    makeplot("rectangle_multiregion")
    makeplot("rectangle_subgrid")
    makeplot("quadrilateral")

    makeplot2("plotting_grid3d")
    makeplot2("plotting_func3d")
    makeplot2("plotting_multiscene")
end


function mkdocs()

    Literate.markdown(examples1d, example_md_dir, documenter=false,info=false)
    Literate.markdown(examples2d, example_md_dir, documenter=false,info=false)
    Literate.markdown(examples3d, example_md_dir, documenter=false,info=false)
    Literate.markdown(plotting, example_md_dir, documenter=false,info=false)
    
    generated_examples=joinpath.("examples",filter(x->endswith(x, ".md"),readdir(example_md_dir)))

    makeplots(example_md_dir)
    
    makedocs(sitename="ExtendableGrids.jl",
    modules = [ExtendableGrids],
    doctest = false,
    clean = true,
             authors = "J. Fuhrmann, Ch. Merdon, J. Krummbiegel",
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
                 "more.md",
                 "plot.md",
                 "tokenstream.md",
                 "allindex.md",
                 "Examples" => generated_examples
             ])
end

mkdocs()

deploydocs(repo = "github.com/j-fu/ExtendableGrids.jl.git")

