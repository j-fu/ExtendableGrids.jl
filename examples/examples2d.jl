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
    rect!(grid,[0.25,0.25],[0.75,0.75];region=2, bregion=5)
    subgrid(grid,[1])
end
# ![](rectangle_subgrid.svg)


# 
# ## Rect2d with bregion function
#
# Here, we use function as bregion parameter - this allows to
# have no bfaces at the interface between the two rects.
function rect2d_bregion_function()
    X=collect(0:0.5:10)
    Y=collect(0:0.5:10)
    grid=simplexgrid(X,Y)
    rect!(grid,[5,4],[9,6];region=2, bregions=[5,5,5,5])

    rect!(grid,[4,2],[5,8];region=2, bregion= cur-> cur == 5  ? 0 : 8   )
    
    subgrid(grid,[2])
    
end
# ![](rect2d_bregion_function.svg)
