# 
# # 301: Region editor and subgrids
# ([source code](SOURCE_URL))
# 

module Example401_RegionEdit_Subgrid

using ExtendableGrids


function main1d(;Plotter=nothing)
    X=collect(0:0.01:1)
    grid=simplexgrid(X)
    bfacemask!(grid,[0.5],[0.5],3)
    cellmask!(grid,[0.0],[0.25],2)
    cellmask!(grid,[0.20],[0.5],3)
    sub=subgrid(grid,[2,3])
    ExtendableGrids.plot(sub,Plotter=Plotter)
    (num_nodes(sub),num_cells(sub),num_bfaces(sub))==(51,50,2)
end

function main2d(;Plotter=nothing)
    X=collect(0:0.05:1)
    Y=collect(0:0.05:1)
    grid=simplexgrid(X,Y)
    cellmask!(grid,[0.25,0.25],[0.75,0.75],2)
    bfacemask!(grid,[0.25,0.25],[0.25,0.75],5)
    bfacemask!(grid,[0.25,0.25],[0.75,0.25],5)
    bfacemask!(grid,[0.25,0.75],[0.75,0.75],5)
    bfacemask!(grid,[0.75,0.25],[0.75,0.75],5)
    sub=subgrid(grid,[1])
    ExtendableGrids.plot(sub,Plotter=Plotter)
    (num_nodes(sub),num_cells(sub),num_bfaces(sub))==(360, 600, 120)
end

function test()
    main1d()&&main2d()
end
end
