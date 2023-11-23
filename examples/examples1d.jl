# 1D Grid examples
# ===============
using ExtendableGrids

# ## Interval from vector
#
function interval_from_vector()
    X = collect(0:0.05:1)
    grid = simplexgrid(X)
end
# ![](interval_from_vector.svg)

# 
# ## Interval with local refinement
# 
function interval_localref()
    XLeft = geomspace(0.0, 0.5, 0.1, 0.01)
    XRight = geomspace(0.5, 1.0, 0.01, 0.1)
    X = glue(XLeft, XRight)
    grid = simplexgrid(X)
end
# ![](interval_localref.svg)

#
# ## Interval with  multiple regions
#
function interval_multiregion()
    X = collect(0:0.05:1)
    grid = simplexgrid(X)
    cellmask!(grid, [0.0], [0.5], 3)
    bfacemask!(grid, [0.5], [0.5], 4)
    grid
end
# ![](interval_multiregion.svg)
#
# ## Multiple regions and subgrid
#
function interval_subgrid()
    X = collect(0:0.01:1)
    grid = simplexgrid(X)
    bfacemask!(grid, [0.5], [0.5], 3)
    cellmask!(grid, [0.0], [0.25], 2)
    cellmask!(grid, [0.20], [0.5], 3)
    subgrid(grid, [2, 3])
end
# ![](interval_subgrid.svg)
# ## CI callbacks
# Unit tests
using Test

function runtests()
    @test numbers_match(interval_from_vector(), 21, 20, 2)
    @test numbers_match(interval_localref(), 27, 26, 2)
    @test numbers_match(interval_multiregion(), 21, 20, 3)
    @test numbers_match(interval_subgrid(), 51, 50, 2)
end

# Plot generation
using GridVisualize
function generateplots(picdir; Plotter = nothing)
    if isdefined(Plotter, :Makie)
        size = (500, 200)
        legend = :rt
        Plotter.save(joinpath(picdir, "interval_from_vector.svg"), gridplot(interval_from_vector(); Plotter, size, legend))
        Plotter.save(joinpath(picdir, "interval_localref.svg"), gridplot(interval_localref(); Plotter, size, legend))
        Plotter.save(joinpath(picdir, "interval_multiregion.svg"), gridplot(interval_multiregion(); Plotter, size, legend))
        Plotter.save(joinpath(picdir, "interval_subgrid.svg"), gridplot(interval_subgrid(); Plotter, size, legend))
    end
end
