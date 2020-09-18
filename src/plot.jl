using .GridRosetta

"""
$(TYPEDSIGNATURES)
Check if Plotter is VTKView
"""
isvtkview(Plotter)= (typeof(Plotter)==Module)&&isdefined(Plotter,:StaticFrame)

"""
$(TYPEDSIGNATURES)
Check if Plotter is PyPlot
"""
ispyplot(Plotter)= (typeof(Plotter)==Module)&&isdefined(Plotter,:Gcf)

"""
$(TYPEDSIGNATURES)
Check if Plotter is Plots
"""
isplots(Plotter)= (typeof(Plotter)==Module) && isdefined(Plotter,:gr)


"""
$(TYPEDSIGNATURES)
Check if Plotter is Makie
"""
ismakie(Plotter)= (typeof(Plotter)==Module)&&isdefined(Plotter,:AbstractPlotting)


"""
$(TYPEDSIGNATURES)
Return tridata to be splatted to PyPlot calls
"""
function tridata(grid)
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    coord[1,:], coord[2,:],transpose(cellnodes.-1)
end

"""
$(TYPEDSIGNATURES)
Return rectangular grid data + function to be splatted into Plots calls
"""
function rectdata(grid,U)
    if dim_grid(grid)==1 && haskey(grid,XCoordinates) 
        return grid[XCoordinates],U
    end
    if dim_grid(grid)==2 && haskey(grid,XCoordinates) && haskey(grid,YCoordinates)
        X=grid[XCoordinates]
        Y=grid[YCoordinates]
        return X,Y,transpose(reshape(U,length(X),length(Y)))
    end
    error("no rectdata on grid")
    nothing
end


"""
$(TYPEDSIGNATURES)
Color scale for grid colors.
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

Plot grid. 

Plotter defaults to `nothing` and can be `PyPLot`, `Plots`, `VTKView`.
"""
function plot(grid::ExtendableGrid;
              Plotter=nothing,
              aspect=1,
              clear=true,
              show=true,
              legend=(1.2,0.5),
              p=nothing)

    
    if ismakie(Plotter)
        mesh=Mesh(grid)
        if p==nothing
            scene=Plotter.Scene()
        elseif typeof(p)==Plotter.Scene
            scene=p
        elseif typeof(p)<:Tuple
            node=p[2]
            node[]=mesh
            return p
        end
        node=Plotter.Node(mesh)
        Plotter.wireframe!(p,Plotter.lift(a->a,node))
        return (scene,node)
    end
    
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

        cellregions=grid[CellRegions]
        cellnodes=grid[CellNodes]
        coord=grid[Coordinates]
        ncellregions=grid[NumCellRegions]
        bfacenodes=grid[BFaceNodes]
        bfaceregions=grid[BFaceRegions]
        nbfaceregions=grid[NumBFaceRegions]
        ncellregions=grid[NumCellRegions]

        crflag=ones(Bool,ncellregions)
        brflag=ones(Bool,nbfaceregions)
        
        if dim_space(grid)==1
            xmin=minimum(coord)
            xmax=maximum(coord)
            h=(xmax-xmin)/40.0
            ax.set_aspect(1)
            ax.get_yaxis().set_ticks([])
            ax.set_ylim(-5*h,xmax-xmin)
            for icell=1:num_cells(grid)
                ireg=cellregions[icell]
                label = crflag[ireg] ? "cellregion $(ireg)" : ""
                crflag[ireg]=false
                
                rgb=frgb(Plotter,ireg,ncellregions)
                x1=coord[1,cellnodes[1,icell]]
                x2=coord[1,cellnodes[2,icell]]
                Plotter.plot([x1,x1],[-h,h],linewidth=0.5,color="k",label="")
                Plotter.plot([x2,x2],[-h,h],linewidth=0.5,color="k",label="")
                Plotter.plot([x1,x2],[0,0],linewidth=3.0,color=rgb,label=label)
            end
            
            for ibface=1:num_bfaces(grid)
                ireg=bfaceregions[ibface]
                if ireg >0
                    label = brflag[ireg] ? "boundary $(ireg)" : ""
                    brflag[ireg]=false
                    rgb=frgb(Plotter,bfaceregions[ibface],nbfaceregions)
                    x1=coord[1,bfacenodes[1,ibface]]
                    Plotter.plot([x1,x1],[-2*h,2*h],linewidth=3.0,color=rgb,label=label)
                end
            end
            Plotter.legend()
        end



        if dim_space(grid)==2
            ax.set_aspect(1)
            tridat=tridata(grid)
            Plotter.tripcolor(tridat...,facecolors=grid[CellRegions],cmap="Pastel2")
            cbar=Plotter.colorbar(ticks=collect(1:ncellregions))
            Plotter.triplot(tridat...,color="k",linewidth=0.5)
            # see https://gist.github.com/gizmaa/7214002
            xc=[coord[:,bfacenodes[1,i]] for i=1:size(bfacenodes,2)]
            yc=[coord[:,bfacenodes[2,i]] for i=1:size(bfacenodes,2)]
            rgb=[frgb(Plotter,bfaceregions[i],nbfaceregions) for i=1:length(bfaceregions)]
            ax.add_collection(Plotter.matplotlib.collections.LineCollection(collect(zip(xc,yc)),colors=rgb,linewidth=3))

            for i=1:nbfaceregions
                Plotter.plot(coord[:,1], coord[:,1],label="boundary $(i)", color=frgb(Plotter,i,nbfaceregions))
            end
            Plotter.legend(loc=legend)
        end
    end

    if isplots(Plotter)
        if p==nothing
            p=Plotter.plot()
        end

        cellregions=grid[CellRegions]
        cellnodes=grid[CellNodes]
        coord=grid[Coordinates]
        ncellregions=grid[NumCellRegions]
        bfacenodes=grid[BFaceNodes]
        bfaceregions=grid[BFaceRegions]
        nbfaceregions=grid[NumBFaceRegions]

        if dim_space(grid)==1
            xmin=minimum(coord)
            xmax=maximum(coord)
            h=(xmax-xmin)/20.0
            
            for icell=1:num_cells(grid)
                rgb=frgb(Plotter,cellregions[icell],ncellregions)
                x1=coord[1,cellnodes[1,icell]]
                x2=coord[1,cellnodes[2,icell]]
                Plotter.plot!(p,[x1,x1],[-h,h],linewidth=0.5,color=:black,label="")
                Plotter.plot!(p,[x2,x2],[-h,h],linewidth=0.5,color=:black,label="")
                Plotter.plot!(p,[x1,x2],[0,0],linewidth=3.0,color=rgb,label="")
            end
            
            for ibface=1:num_bfaces(grid)
                if bfaceregions[ibface]>0
                rgb=frgb(Plotter,bfaceregions[ibface],nbfaceregions)
                    x1=coord[1,bfacenodes[1,ibface]]
                    Plotter.plot!(p,[x1,x1],[-2*h,2*h],linewidth=3.0,color=rgb,label="")
                end
            end
        end
        
        if dim_space(grid)==2
            for icell=1:num_cells(grid)
                rgb=frgb(Plotter,cellregions[icell],ncellregions,pastel=true)
                inode1=cellnodes[1,icell]
                inode2=cellnodes[2,icell]
                inode3=cellnodes[3,icell]
                # https://github.com/JuliaPlotter/Plotter.jl/issues/605
                
                tri=Plotter.Shape([coord[1,inode1],coord[1,inode2], coord[1,inode3]],[coord[2,inode1],coord[2,inode2],coord[2,inode3]])
                Plotter.plot!(p,tri,color=rgb,label="")
            end
            for icell=1:num_cells(grid)
                inode1=cellnodes[1,icell]
                inode2=cellnodes[2,icell]
                inode3=cellnodes[3,icell]
                Plotter.plot!(p, [coord[1,inode1],coord[1,inode2]],[coord[2,inode1],coord[2,inode2]]  ,linewidth=0.5,color=:black,label="")
                Plotter.plot!(p, [coord[1,inode1],coord[1,inode3]],[coord[2,inode1],coord[2,inode3]]  ,linewidth=0.5,color=:black,label="")
                Plotter.plot!(p, [coord[1,inode2],coord[1,inode3]],[coord[2,inode2],coord[2,inode3]]  ,linewidth=0.5,color=:black,label="")
            end
            for ibface=1:num_bfaces(grid)
                rgb=frgb(Plotter,bfaceregions[ibface],nbfaceregions)
                inode1=bfacenodes[1,ibface]
                inode2=bfacenodes[2,ibface]
                Plotter.plot!(p,[coord[1,inode1],coord[1,inode2]],[coord[2,inode1],coord[2,inode2]]  ,linewidth=5,color=rgb,label="")
            end
        end
        if show
            Plotter.gui(p)
        end
        return p
    end


end

"""
$(TYPEDSIGNATURES)

Plot lowest order continuous finite element function on grid
defined by values on the grid nodes.


Keyword arguments:
- Plotter: defaults to `nothing` and can be `PyPLot`, `Plots`, `VTKView`.
- color:  color of plot on 1D grid
- cmap:  color map for heatmap plot
- label: label of plot
- levels: number of isolevels
- aspect: aspect ratio of plot
- clear: if true (default) clear plot before plotting
- show: if true (default) show plot
"""
function plot(grid::ExtendableGrid, U::AbstractVector;
              Plotter=nothing,
              color=(0,0,0),
              cmap="hot",
              label="",
              colorlevels=51,
              isolines=11,
              aspect=1,
              clear=true,
              show=true,
              cbar=true,
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
        PyPlot=Plotter
        if clear
            PyPlot.clf()
        end
        if dim_space(grid)==1
            for icell=1:num_cells(grid)
                i1=cellnodes[1,icell]
                i2=cellnodes[2,icell]
                x1=coord[1,i1]
                x2=coord[1,i2]
                if icell==1 && label !=""
                    PyPlot.plot([x1,x2],[U[i1],U[i2]],color=color,label=label)
                else
                    PyPlot.plot([x1,x2],[U[i1],U[i2]],color=color)
                end                
            end
        end
        
        if dim_space(grid)==2
            ax=PyPlot.matplotlib.pyplot.gca()
            ax.set_aspect(aspect)
            umin=minimum(U)
            umax=maximum(U)
            if typeof(colorlevels)<:Number
                colorlevels=collect(umin:(umax-umin)/(colorlevels-1):umax)
            end
            if typeof(isolines)<:Number
                isolines=collect(umin:(umax-umin)/(isolines-1):umax)
            end
            PyPlot.tricontourf(tridata(grid)...,U;levels=colorlevels,cmap=PyPlot.ColorMap(cmap))
            if cbar
                PyPlot.colorbar(ticks=isolines,boundaries=colorlevels)
            end
            PyPlot.tricontour(tridata(grid)...,U,colors="k",levels=isolines)
            if show
                PyPlot.pause(1.0e-10)
                PyPlot.show()
            end
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

function plot(gf::GridFactory; Plotter=nothing,title="")
    
    if ispyplot(Plotter)
        triin=nothing
        try
            triin=triangulateio(gf)
        catch err
            @error "Incomplete geometry description"
            rethrow(err)
        end
        if typeof(gf.unsuitable)!=Nothing
            triunsuitable(gf.unsuitable)
        end
        triout,vorout=Triangulate.triangulate(gf.flags,triin)
        PyPlot=Plotter
        PyPlot.clf()
        PyPlot.suptitle(title)
        PyPlot.subplot(121)
        PyPlot.title("In")
        Triangulate.plot(PyPlot,triin)
        PyPlot.subplot(122)
        PyPlot.title("Out")
        Triangulate.plot(PyPlot,triout)
        PyPlot.tight_layout()
    end
end


