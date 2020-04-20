isvtkview(Plotter)= (typeof(Plotter)==Module)&&isdefined(Plotter,:StaticFrame)
ispyplot(Plotter)= (typeof(Plotter)==Module)&&isdefined(Plotter,:Gcf)
isplots(Plotter)= (typeof(Plotter)==Module) && isdefined(Plotter,:gr)

function tridata(grid)
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    coord[1,:], coord[2,:],transpose(cellnodes.-1)
end


"""
Plot color scale for grid colors.
"""
function frgb(Plotter,i,max;pastel=false)
    x=Float64(i-1)/Float64(max)
    if (x<0.5)
        r=1.0-2.0*x
        g=2.0*x
        b=0.0
    else
        r=0.0
        g=2.0-2.0*x
        b=2.0*x-1.0
    end
    if pastel
        r=0.5+0.5*r
        g=0.5+0.5*g
        b=0.5+0.5*b
    end
    if ispyplot(Plotter)
        return (r,g,b)
    end
    if isplots(Plotter)
        return Plotter.RGB(r,g,b)
    end
end


function plot(grid::ExtendableGrid; Plotter=nothing)

    
    if isvtkview(Plotter)
        frame=Plotter.StaticFrame()
        Plotter.clear!(frame)
        dataset=Plotter.DataSet()
        Plotter.simplexgrid!(dataset,grid[Coordinates],grid[CellNodes])
        Plotter.boundarygrid!(dataset,grid[BFaceNodes])
        Plotter.boundarymarker!(dataset,grid[BFaceRegions])
        Plotter.cellmarker!(dataset,grid[CellRegions])
        gridview=Plotter.GridView()
        Plotter.data!(gridview,dataset)
        Plotter.addview!(frame,gridview)
        Plotter.display(frame)
    end

    if ispyplot(Plotter)
        Plotter.clf()
        ax=Plotter.matplotlib.pyplot.gca()
        ax.set_aspect(1)
        tridat=tridata(grid)
        cellregions=grid[CellRegions]
        coord=grid[Coordinates]
        bfacenodes=grid[BFaceNodes]
        bfaceregions=grid[BFaceRegions]
        Plotter.tripcolor(tridat...,facecolors=grid[CellRegions],cmap="Pastel2")
        Plotter.triplot(tridat...,color="k",linewidth=0.5)
        
        # see https://gist.github.com/gizmaa/7214002
        xc=[coord[:,bfacenodes[1,i]] for i=1:size(bfacenodes,2)]
        yc=[coord[:,bfacenodes[2,i]] for i=1:size(bfacenodes,2)]
        rgb=[frgb(Plotter,bfaceregions[i],maximum(bfaceregions)) for i=1:length(bfaceregions)]
        ax.add_collection(Plotter.matplotlib.collections.LineCollection(collect(zip(xc,yc)),colors=rgb,linewidth=3))
    end
end

