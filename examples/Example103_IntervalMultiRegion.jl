#= 

# 103: Interval with multiple regions
([source code](SOURCE_URL))

=#

module Example103_IntervalMultiRegion
using ExtendableGrids

function main(;plotter=nothing)

    X=collect(0:0.05:1)
    grid=simplexgrid(X)
    cellmask!(grid,[0.0],[0.5],3)
    bfacemask!(grid,[0.5], [0.5],4)
    ExtendableGrids.plot(grid,Plotter=plotter)
    (num_nodes(grid),num_cells(grid),num_bfaces(grid))
end
function test()
    main()==(21,20,3)
end
end
