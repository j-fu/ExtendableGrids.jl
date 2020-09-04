#= 

# 252: Locally refined square using GridFactory
([source code](SOURCE_URL))

=#

module Example252_GridFactorySquareLocalRef
using ExtendableGrids

function main(;plotter=nothing)
    
    factory=GridFactory(dim_space=2)
    appendflags!(factory,"u")
    p1=point!(factory,0,0)
    p2=point!(factory,1,0)
    p3=point!(factory,1,1)
    p4=point!(factory,0,1)
    facet!(factory,p1,p2,region=1)
    facet!(factory,p2,p3,region=2)
    facet!(factory,p3,p4,region=3)
    facet!(factory,p4,p1,region=4)
    
    cellregion!(factory,0.5,0.5,region=1,volume=0.01)

    function unsuitable(x1,y1,x2,y2,x3,y3, area)
        bary_x=x1+x2+x3
        bary_y=y2+y2+y3
        if area > 0.001*bary_x
            return 1
        else
            return 0
        end
    end
    unsuitable!(factory, unsuitable)

    grid=simplexgrid(factory)
    @show grid

    
    ExtendableGrids.plot(grid,Plotter=plotter)
    (num_nodes(grid),num_cells(grid),num_bfaces(grid))
end

function test()
    main()==(2359, 4128, 588)
end
end
