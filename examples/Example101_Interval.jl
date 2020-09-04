#= 

# 101: Interval
([source code](SOURCE_URL))

=#

module Example101_Interval
using ExtendableGrids

function main(;plotter=nothing)

    X=collect(0:0.05:1)
    grid=simplexgrid(X)
    ExtendableGrids.plot(grid,Plotter=plotter)
    (num_nodes(grid),num_cells(grid),num_bfaces(grid))
    
end
function test()
    main()==(21,20,2)
end
end
