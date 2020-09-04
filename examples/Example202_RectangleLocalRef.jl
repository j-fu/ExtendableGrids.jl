#= 

# 202: Rectangle with local refinement
([source code](SOURCE_URL))

=#

module Example202_RectangleLocalRef
using ExtendableGrids

function main(;plotter=Nothing)

    XLeft=geomspace(0.0,0.5,0.1, 0.01)
    XRight=geomspace(0.5,1.0,0.01,0.1)
    X=glue(XLeft, XRight)
    grid=simplexgrid(X,X)
    ExtendableGrids.plot(grid,Plotter=plotter)
    (num_nodes(grid),num_cells(grid),num_bfaces(grid))
    
end
function test()
    main()==(729, 1352, 104)
end
end
