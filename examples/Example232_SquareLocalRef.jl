#= 

# 232: Locally refined square
([source code](SOURCE_URL))

=#

module Example232_SquareLocalRef
using ExtendableGrids

function main(;plotter=Nothing)

    function unsuitable(x1,y1,x2,y2,x3,y3, area)
        bary_x=x1+x2+x3
        bary_y=y2+y2+y3
        if area > 0.001*bary_x
            return 1
        else
            return 0
        end
    end
    
    grid=simplexgrid(points=[0 0 ; 0 1 ; 1 1 ; 1 0]',
                     bfaces=[1 2 ; 2 3 ; 3 4 ; 4 1 ]',
                     bfaceregions=[1, 2, 3, 4],
                     regionpoints=[0.5 0.5;]',
                     regionnumbers=[1],
                     regionvolumes=[0.01],
                     flags="pAaqQDu",
                     unsuitable=unsuitable)
    ExtendableGrids.plot(grid,Plotter=plotter)
    (num_nodes(grid),num_cells(grid),num_bfaces(grid))
    
end
function test()
    main()==(2359, 4128, 588)
end
end
