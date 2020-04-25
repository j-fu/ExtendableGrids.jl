"""
$(SIGNATURES)
Check if Plotter is VTKView
"""
isvtkview(Plotter)= (typeof(Plotter)==Module)&&isdefined(Plotter,:StaticFrame)

"""
$(SIGNATURES)
Check if Plotter is PyPlot
"""
ispyplot(Plotter)= (typeof(Plotter)==Module)&&isdefined(Plotter,:Gcf)

"""
$(SIGNATURES)
Check if Plotter is Plots
"""
isplots(Plotter)= (typeof(Plotter)==Module) && isdefined(Plotter,:gr)

"""
$(SIGNATURES)
Return tridata to be splatted to PyPlot calls
"""
function tridata(grid)
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    coord[1,:], coord[2,:],transpose(cellnodes.-1)
end


"""
$(TYPEDSIGNATURES)
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

"""
$(TYPEDSIGNATURES)

Plot grid
"""
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
        nbfaceregions=grid[NumBFaceRegions]
        Plotter.tripcolor(tridat...,facecolors=grid[CellRegions],cmap="Pastel2")
        Plotter.triplot(tridat...,color="k",linewidth=0.5)
        
        # see https://gist.github.com/gizmaa/7214002
        xc=[coord[:,bfacenodes[1,i]] for i=1:size(bfacenodes,2)]
        yc=[coord[:,bfacenodes[2,i]] for i=1:size(bfacenodes,2)]
        rgb=[frgb(Plotter,bfaceregions[i],nbfaceregions) for i=1:length(bfaceregions)]
        ax.add_collection(Plotter.matplotlib.collections.LineCollection(collect(zip(xc,yc)),colors=rgb,linewidth=3))
    end
end

"""
$(TYPEDSIGNATURES)

Plot lowest order continuous finite element function on grid
"""
function plot(grid::ExtendableGrid, U::AbstractArray;
              Plotter=nothing,
              color=(0,0,0),
              cmap="hot",
              label="",
              levels=10,
              aspect=1,
              clear=true,
              show=true,
              p=nothing)

    cellnodes=grid[CellNodes]
    coord=grid[Coordinates]

    if isvtkview(Plotter)
        if dim_space(grid)==2
            frame=Plotter.StaticFrame()
            Plotter.clear!(frame)
            dataset=Plotter.DataSet()
            Plotter.simplexgrid!(dataset,grid[Coordinates],grid[CellNodes])
            Plotter.pointscalar!(dataset,U,"V")
            scalarview=Plotter.ScalarView()
            Plotter.data!(scalarview,dataset,"V")
            Plotter.addview!(frame,scalarview)
            Plotter.display(frame)
        end        
    end
    
    if ispyplot(Plotter)
        if clear
            Plotter.clf()
        end
        if dim_space(grid)==1
            for icell=1:num_cells(grid)
                i1=cellnodes[1,icell]
                i2=cellnodes[2,icell]
                x1=coord[1,i1]
                x2=coord[1,i2]
                if icell==1 && label !=""
                    Plotter.plot([x1,x2],[U[i1],U[i2]],color=color,label=label)
                else
                    Plotter.plot([x1,x2],[U[i1],U[i2]],color=color)
                end                
            end
        end
        
        if dim_space(grid)==2
            ax=Plotter.matplotlib.pyplot.gca()
            ax.set_aspect(aspect)
            plotted=Plotter.tricontourf(tridata(grid)...,U;levels=levels,cmap=cmap)
            Plotter.tricontour(tridata(grid)...,U,colors="k",levels=levels)
            return plotted
        end
    end

    if isplots(Plotter)
        if p==nothing
            p=Plotter.plot()
        end
        if dim_space(grid)==1
            for icell=1:num_cells(grid)
                i1=cellnodes[1,icell]
                i2=cellnodes[2,icell]
                x1=coord[1,i1]
                x2=coord[1,i2]
                if icell==1 && label !=""
                    Plotter.plot!(p,[x1,x2],[U[i1],U[i2]],linecolor=Plotter.RGB(color...),label=label)
                else
                    Plotter.plot!(p,[x1,x2],[U[i1],U[i2]],linecolor=Plotter.RGB(color...),label="")
                end                
            end
        end
        if dim_space(grid)==2
            println("Not available for Plots, see e.g. https://github.com/JuliaPlots/Plots.jl/issues/392")
        end
        if show
            Plotter.gui(p)
        end
        return p
    end
end
