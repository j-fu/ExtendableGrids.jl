# 3D Grid examples
# ===============
# 

# ## Quadrilateral
function quadrilateral(;hx=0.25, hy=0.2, hz=0.1)
    X=collect(0:hx:1)
    Y=collect(0:hy:1)
    Z=collect(0:hz:1)
    simplexgrid(X,Y,Z)
end
# ![](quadrilateral.svg)


function mask_bedges()
    grid    = quadrilateral(hx=0.25, hy=0.25, hz=0.25)

    bedgemask!(grid, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1)
    bedgemask!(grid, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2)
    bedgemask!(grid, [0.0, 1.0, 0.0], [0.0, 1.0, 1.0], 3)
    bedgemask!(grid, [0.0, 0.0, 1.0], [0.0, 1.0, 1.0], 4)
    bedgemask!(grid, [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 5)

    true
end


