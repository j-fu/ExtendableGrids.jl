# Gmsh examples
# ===============

using ExtendableGrids
using Gmsh: gmsh

function gmsh2d()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("t1")

    lc = 1e-2
    gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
    gmsh.model.geo.addPoint(0.1, 0, 0, lc, 2)
    gmsh.model.geo.addPoint(0.1, 0.3, 0, lc, 3)

    p4 = gmsh.model.geo.addPoint(0, 0.3, 0, lc)

    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(3, 2, 2)
    gmsh.model.geo.addLine(3, p4, 3)
    gmsh.model.geo.addLine(4, 1, p4)

    gmsh.model.geo.addCurveLoop([4, 1, -2, 3], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)

    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(0, [1, 2], 1)
    gmsh.model.addPhysicalGroup(1, [1, 2], 2)
    gmsh.model.addPhysicalGroup(2, [1], 6)

    gmsh.model.setPhysicalName(2, 6, "My surface")

    gmsh.model.mesh.generate(2)
    grid = ExtendableGrids.simplexgrid_from_gmsh(gmsh.model)
    gmsh.finalize()
    grid
end
# ![](gmsh2d.svg)

# ## CI callbacks for [ExampleJuggler.jl](https://github.com/j-fu/ExampleJuggler.jl)
# Unit tests
using Test
function runtests()
    grid = gmsh2d()
    @test num_nodes(grid) > 0 && num_cells(grid) > 0 && num_bfaces(grid) > 0
end

# Plot generation
using GridVisualize
function generateplots(picdir; Plotter = nothing)
    if isdefined(Plotter, :Makie)
        size = (300, 300)
        Plotter.save(joinpath(picdir, "gmsh2d.svg"), gridplot(gmsh2d(); Plotter, size))
    end
end
