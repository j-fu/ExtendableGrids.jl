#= 

# 233: Square with hole
([source code](SOURCE_URL))

=#

module Example233_SquareWithHole
using ExtendableGrids

function main(;plotter=nothing)

    function unsuitable(x1,y1,x2,y2,x3,y3, area)
        bary_x=x1+x2+x3
        bary_y=y2+y2+y3
        if area > 0.001*bary_x
            return 1
        else
            return 0
        end
    end
    
    grid=simplexgrid(points=[0 0 ; 0 1 ; 1 1 ; 1 0; 0.3 0.3 ; 0.3 0.7 ; 0.7 0.7  ; 0.7 0.3]',
                     bfaces=[1 2 ; 2 3 ; 3 4 ; 4 1 ; 5 6 ; 6 7; 7 8; 8 5]',
                     bfaceregions=[1, 2, 3, 4,5,5,5,5],
                     regionpoints=[0.25 0.25;0.5 0.5]',
                     regionnumbers=[1, 0],
                     regionvolumes=[0.01, 1],
                     flags="pAaqQD")

    ExtendableGrids.plot(grid,Plotter=plotter)
    (num_nodes(grid),num_cells(grid),num_bfaces(grid))
end
function test()
    main()==(93, 139, 47)
end
end
