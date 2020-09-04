#= 

# 253: Square with hole via GridFactory
([source code](SOURCE_URL))

=#

module Example253_GridFactorySquareWithHole
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


    h1=point!(factory,0.3, 0.3)
    h2=point!(factory,0.3, 0.7)
    h3=point!(factory,0.7, 0.7)
    h4=point!(factory,0.7, 0.3)

    facet!(factory,h1,h2,region=5)
    facet!(factory,h2,h3,region=5)
    facet!(factory,h3,h4,region=5)
    facet!(factory,h4,h1,region=5)

    hole!(factory, 0.5, 0.5)
    cellregion!(factory,0.25,0.25,region=1,volume=0.01)

    grid=simplexgrid(factory)
    
    ExtendableGrids.plot(grid,Plotter=plotter)
    (num_nodes(grid),num_cells(grid),num_bfaces(grid))
end

function test()
    main()==(93, 139, 47)
end
end
