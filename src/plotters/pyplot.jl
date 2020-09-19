function initialize_context!(ctx::PlotterContext,::Type{PyPlotType})
    PyPlot=ctx[:Plotter]
    ctx
end

function prepare_figure!(ctx::PlotterContext)
    PyPlot=ctx[:Plotter]
    if !haskey(ctx,:figure)
        res=ctx[:resolution]
        ctx[:figure]=PyPlot.figure(ctx[:fignumber],figsize=(res[1]/100,res[2]/100),dpi=100)
    end
    if ctx[:clear]
        PyPlot.clf()
    end
    ctx
end


"""
$(TYPEDSIGNATURES)
Return tridata to be splatted to PyPlot calls
"""
function tridata(grid)
    coord=grid[Coordinates]
    cellnodes=Matrix(grid[CellNodes])
    coord[1,:], coord[2,:],transpose(cellnodes.-1)
end


function plot!(ctx, ::Type{PyPlotType}, ::Type{Val{1}}, grid)
    PyPlot=ctx[:Plotter]

    prepare_figure!(ctx)

    ax=PyPlot.matplotlib.pyplot.gca()
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
        
        rgb=frgb(PyPlot,ireg,ncellregions)
        x1=coord[1,cellnodes[1,icell]]
        x2=coord[1,cellnodes[2,icell]]
        PyPlot.plot([x1,x1],[-h,h],linewidth=0.5,color="k",label="")
        PyPlot.plot([x2,x2],[-h,h],linewidth=0.5,color="k",label="")
        PyPlot.plot([x1,x2],[0,0],linewidth=3.0,color=rgb,label=label)
    end
    
    for ibface=1:num_bfaces(grid)
        ireg=bfaceregions[ibface]
        if ireg >0
            label = brflag[ireg] ? "boundary $(ireg)" : ""
            brflag[ireg]=false
            rgb=frgb(PyPlot,bfaceregions[ibface],nbfaceregions)
            x1=coord[1,bfacenodes[1,ibface]]
            PyPlot.plot([x1,x1],[-2*h,2*h],linewidth=3.0,color=rgb,label=label)
        end
    end
    PyPlot.legend()
    ctx[:figure]
end



function plot!(ctx, ::Type{PyPlotType}, ::Type{Val{2}},grid)
    PyPlot=ctx[:Plotter]

    prepare_figure!(ctx)

    ax=PyPlot.matplotlib.pyplot.gca()
    cellregions=grid[CellRegions]
    cellnodes=grid[CellNodes]
    coord=grid[Coordinates]
    ncellregions=grid[NumCellRegions]
    nbfaceregions=grid[NumBFaceRegions]
    ncellregions=grid[NumCellRegions]
    if nbfaceregions>0
        bfacenodes=grid[BFaceNodes]
        bfaceregions=grid[BFaceRegions]
    end
    
    crflag=ones(Bool,ncellregions)
    brflag=ones(Bool,nbfaceregions)
    ax.set_aspect(1)
    tridat=tridata(grid)
    PyPlot.tripcolor(tridat...,facecolors=grid[CellRegions],cmap="Pastel2")
    cbar=PyPlot.colorbar(ticks=collect(1:ncellregions))
    PyPlot.triplot(tridat...,color="k",linewidth=0.5)



    if nbfaceregions>0
        # see https://gist.github.com/gizmaa/7214002
        xc=[coord[:,bfacenodes[1,i]] for i=1:size(bfacenodes,2)]
        yc=[coord[:,bfacenodes[2,i]] for i=1:size(bfacenodes,2)]
        rgb=[frgb(PyPlot,bfaceregions[i],nbfaceregions) for i=1:length(bfaceregions)]
        ax.add_collection(PyPlot.matplotlib.collections.LineCollection(collect(zip(xc,yc)),colors=rgb,linewidth=3))
    
        for i=1:nbfaceregions
            PyPlot.plot(coord[:,1], coord[:,1],label="boundary $(i)", color=frgb(PyPlot,i,nbfaceregions))
        end
    end
    PyPlot.legend(loc=ctx[:legend_location])
    PyPlot.show()
    PyPlot.pause(0.0001)
    ctx[:figure]
end


function plot!(ctx, ::Type{PyPlotType}, ::Type{Val{1}},grid, func)
    PyPlot=ctx[:Plotter]

    prepare_figure!(ctx)

    cellnodes=grid[CellNodes]
    coord=grid[Coordinates]
    for icell=1:num_cells(grid)
        i1=cellnodes[1,icell]
        i2=cellnodes[2,icell]
        x1=coord[1,i1]
        x2=coord[1,i2]
        if icell==1 && ctx[:label] !=""
            PyPlot.plot([x1,x2],[func[i1],func[i2]],color=ctx[:color],label=ctx[:label])
        else
            PyPlot.plot([x1,x2],[func[i1],func[i2]],color=ctx[:color])
        end                
    end
    if ctx[:legend]
        PyPlot.legend(loc=ctx[:legend_location])
    end
    if ctx[:axisgrid]
        PyPlot.grid()
    end
    ctx[:figure]
end

function plot!(ctx, ::Type{PyPlotType}, ::Type{Val{2}},grid, func)
    PyPlot=ctx[:Plotter]

    prepare_figure!(ctx)

    ax=PyPlot.matplotlib.pyplot.gca()
    ax.set_aspect(ctx[:aspect])
    umin=minimum(func)
    umax=maximum(func)
    colorlevels=collect(umin:(umax-umin)/(ctx[:colorlevels]-1):umax)
    isolines=collect(umin:(umax-umin)/(ctx[:isolines]-1):umax)
    if !haskey(ctx,:grid) || !seemingly_equal(ctx[:grid],grid)
        ctx[:grid]=grid
        ctx[:tridata]=tridata(grid)
    end
    PyPlot.tricontourf(ctx[:tridata]...,func;levels=colorlevels,cmap=PyPlot.ColorMap(ctx[:colormap]))
    if ctx[:colorbar]
        PyPlot.colorbar(ticks=isolines,boundaries=colorlevels)
    end
    PyPlot.tricontour(ctx[:tridata]...,func,colors="k",levels=isolines)
    PyPlot.show()
    PyPlot.pause(0.0001)
    ctx[:figure]
end


function plot!(ctx, ::Type{PyPlotType}, gf::GridFactory)
    PyPlot=ctx[:Plotter]

    prepare_figure!(ctx)

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
    PyPlot.subplot(121)
    PyPlot.title("In")
    Triangulate.plot(PyPlot,triin)
    PyPlot.subplot(122)
    PyPlot.title("Out")
    Triangulate.plot(PyPlot,triout)
    PyPlot.tight_layout()
    ctx[:figure]
end

