#= 

# 231: Square
([source code](SOURCE_URL))

=#

module Example231_Square
using ExtendableGrids

function main(;plotter=Nothing)

    grid=simplexgrid(points=[0 0 ; 0 1 ; 1 1 ; 1 0]',
                     bfaces=[1 2 ; 2 3 ; 3 4 ; 4 1 ]',
                     bfaceregions=[1, 2, 3, 4],
                     regionpoints=[0.5 0.5;]',
                     regionnumbers=[1],
                     regionvolumes=[0.01],
                     flags="pAaqQD")
    ExtendableGrids.plot(grid,Plotter=plotter)
    (num_nodes(grid),num_cells(grid),num_bfaces(grid))
    
end
function test()
    main()==(89,144,32)
end
end
