#= 

# 251: Square using GridFactory
([source code](SOURCE_URL))

=#

module Example251_GridFactorySquare
using ExtendableGrids

function main(;plotter=nothing)
    
    factory=GridFactory(dim_space=2)
    
    p1=point!(factory,0,0)
    p2=point!(factory,1,0)
    p3=point!(factory,1,1)
    p4=point!(factory,0,1)
    facet!(factory,p1,p2,region=1)
    facet!(factory,p2,p3,region=2)
    facet!(factory,p3,p4,region=3)
    facet!(factory,p4,p1,region=4)
    
    cellregion!(factory,0.5,0.5,region=1,volume=0.01)
    grid=simplexgrid(factory)
    @show grid
    
    ExtendableGrids.plot(grid,Plotter=plotter)
    (num_nodes(grid),num_cells(grid),num_bfaces(grid))
end

function test()
    main()==(89,144,32)
end
end
