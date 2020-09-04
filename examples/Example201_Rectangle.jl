#= 

# 201: Rectangle
([source code](SOURCE_URL))

=#

module Example201_Rectangle
using ExtendableGrids

function main(;plotter=Nothing)

    X=collect(0:0.05:1)
    Y=collect(0:0.05:1)
    grid=simplexgrid(X,Y)
    ExtendableGrids.plot(grid,Plotter=plotter,legend=(1.1,0.75))
    (num_nodes(grid),num_cells(grid),num_bfaces(grid))
end
function test()
    main()==(441,800,80)
end
end
