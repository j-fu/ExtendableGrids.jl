# 3D Grid examples
# ===============
# 

# ## Quadrilateral
function quadrilateral()
    X=collect(0:0.25:1)
    Y=collect(0:0.2:1)
    Z=collect(0:0.1:1)
    simplexgrid(X,Y,Z)
end
# ![](quadrilateral.svg)


