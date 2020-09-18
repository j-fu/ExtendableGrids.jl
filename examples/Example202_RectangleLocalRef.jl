#= 

# 202: Rectangle with local refinement
([source code](SOURCE_URL))

=#

module Example202_RectangleLocalRef
using ExtendableGrids

function main(;hmin=0.01, hmax=0.1, plotter=nothing, scene=nothing)

    XLeft=geomspace(0.0,0.5,hmax,hmin)
    XRight=geomspace(0.5,1.0,hmin,hmax)
    X=glue(XLeft, XRight)
    grid=simplexgrid(X,X)
    if plotter==nothing
        return (num_nodes(grid),num_cells(grid),num_bfaces(grid))
    end
    ExtendableGrids.plot(grid,Plotter=plotter, p=scene)
    
end
function test()
    main()==(729, 1352, 104)
end
end
