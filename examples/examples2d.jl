# 2D Grid examples
# ===============
# 

# ## Rectangle
function rectangle()
    X=collect(0:0.05:1)
    Y=collect(0:0.05:1)
    simplexgrid(X,X)
end
# ![](rectangle.svg)
# 
# ## Rectangle with local refinement
# 
function rectangle_localref()
    hmin=0.01
    hmax=0.1
    XLeft=geomspace(0.0,0.5,hmax,hmin)
    XRight=geomspace(0.5,1.0,hmin,hmax)
    X=glue(XLeft, XRight)
    simplexgrid(X,X)
end
# ![](rectangle_localref.svg)


# 
# ## Rectangle with multiple regions
# 
function rectangle_multiregion()
    X=collect(0:0.05:1)
    Y=collect(0:0.05:1)
    grid=simplexgrid(X,Y)
    cellmask!(grid,[0.0,0.0],[1.0,0.5],3)
    bfacemask!(grid,[0.0,0.0],[0.0,0.5],5)
    bfacemask!(grid,[1.0,0.0],[1.0,0.5],6)
    bfacemask!(grid,[0.0,0.5],[1.0,0.5],7)
end
# ![](rectangle_multiregion.svg)



# 
# ## Subgrid from rectangle
# 
function rectangle_subgrid()
    X=collect(0:0.05:1)
    Y=collect(0:0.05:1)
    grid=simplexgrid(X,Y)
    cellmask!(grid,[0.25,0.25],[0.75,0.75],2)
    bfacemask!(grid,[0.25,0.25],[0.25,0.75],5)
    bfacemask!(grid,[0.25,0.25],[0.75,0.25],5)
    bfacemask!(grid,[0.25,0.75],[0.75,0.75],5)
    bfacemask!(grid,[0.75,0.25],[0.75,0.75],5)
    subgrid(grid,[1])
end
# ![](rectangle_subgrid.svg)
