# 
# # 203: Rectangle with multiple regions
# ([source code](SOURCE_URL))
# 


module Example203_RectangleMultiRegion

using ExtendableGrids
  
function main(;plotter=nothing)

    X=collect(0:0.05:1)
    Y=collect(0:0.05:1)
    grid=simplexgrid(X,Y)
    cellmask!(grid,[0.0,0.0],[1.0,0.5],3)
    bfacemask!(grid,[0.0,0.0],[0.0,0.5],5)
    bfacemask!(grid,[1.0,0.0],[1.0,0.5],6)
    ExtendableGrids.plot(grid,Plotter=plotter)
    (num_nodes(grid),num_cells(grid),num_bfaces(grid))
end
function test()
    main()==(441,800,80)
end
end
