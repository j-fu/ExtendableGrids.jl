function initialize_gridplot!(p::GridPlotContext,::Type{VTKViewType})
    pctx=p.context
    VTKView=pctx[:Plotter]
    frame=VTKView.StaticFrame()
    pctx[:frame]=frame
    if pctx[:clear]
        VTKView.clear!(frame)
    end
    layout=pctx[:layout]
    tlayout=(layout[2],layout[1])
    for I in CartesianIndices(layout)
        ctx=p.subplots[I]
        ctx[:frame]=frame
    end
    VTKView.layout!(frame,tlayout...)
    VTKView.size!(frame,pctx[:resolution]...)
    pctx
end


function save(fname,p,::Type{VTKViewType})
    VTKView=p.context[:Plotter]
    base,ext=splitext(fname)
    if ext!=".png"
        error("VTKView can only save png files")
    end
    VTKView.writepng(p.context[:frame],fname)
end



function reveal(p::GridPlotContext,::Type{VTKViewType})
    VTKView=p.Plotter
    VTKView.display(p.context[:frame])
end



function reveal(ctx::SubPlotContext,TP::Type{VTKViewType})
    if ctx[:show]||ctx[:reveal]
        reveal(ctx[:GridPlotContext],TP)
    end
end




function gridplot!(ctx, TP::Type{VTKViewType},grid)
    VTKView=ctx[:Plotter]
    frame=ctx[:frame]
    if !haskey(ctx,:dataset)
        dataset=VTKView.DataSet()
        ctx[:dataset]=dataset
    end
    VTKView.simplexgrid!(ctx[:dataset],grid[Coordinates],grid[CellNodes])
    VTKView.boundarygrid!(ctx[:dataset],grid[BFaceNodes])
    VTKView.boundarymarker!(ctx[:dataset],grid[BFaceRegions])
    VTKView.cellmarker!(ctx[:dataset],grid[CellRegions])
    if !haskey(ctx,:gridview)
        ctx[:gridview]=VTKView.GridView()
        VTKView.data!(ctx[:gridview],ctx[:dataset])
        VTKView.addview!(frame,ctx[:gridview],ctx[:iplot]...)
    end
    reveal(ctx,TP)
end


gridplot!(ctx, T::Type{VTKViewType}, ::Type{Val{2}},grid)=gridplot!(ctx, T,grid)
gridplot!(ctx, T::Type{VTKViewType}, ::Type{Val{3}},grid)=gridplot!(ctx, T,grid)


function gridplot!(ctx, TP::Type{VTKViewType},grid,func)
    VTKView=ctx[:Plotter]
    frame=ctx[:frame]
    if !haskey(ctx,:dataset)
        ctx[:dataset]=VTKView.DataSet()
    end
    if !haskey(ctx,:grid)
        ctx[:grid]=grid
        VTKView.simplexgrid!(ctx[:dataset],grid[Coordinates],grid[CellNodes])
    end
    if !haskey(ctx,:scalarview)
        ctx[:scalarview]=VTKView.ScalarView()
        VTKView.addview!(frame,ctx[:scalarview],ctx[:iplot])
        VTKView.data!(ctx[:scalarview],ctx[:dataset],ctx[:label])
    end
    if !seemingly_equal(grid,ctx[:grid])
        VTKView.simplexgrid!(ctx[:dataset],grid[Coordinates],grid[CellNodes])
    end
    VTKView.pointscalar!(ctx[:dataset],func,ctx[:label])
    if dim_space(grid)==3
        VTKView.isolevels!(ctx[:scalarview],[ctx[:flevel]])
        VTKView.show_isosurfaces!(ctx[:scalarview],true)
    end
    reveal(ctx,TP)
end



gridplot!(ctx, T::Type{VTKViewType}, ::Type{Val{2}},grid, func)=gridplot!(ctx, T,grid,func)
gridplot!(ctx, T::Type{VTKViewType}, ::Type{Val{3}},grid, func)=gridplot!(ctx, T,grid,func)
gridplot!(ctx, T::Type{VTKViewType}, ::Type{Val{1}},grid)=nothing

function gridplot!(ctx, TP::Type{VTKViewType}, ::Type{Val{1}},grid, func)
    VTKView=ctx[:Plotter]
    frame=ctx[:frame]
    if !haskey(ctx,:plot)
        ctx[:plot]=VTKView.XYPlot()
        VTKView.addview!(frame,ctx[:plot],ctx[:iplot])
        VTKView.xrange!(ctx[:plot],ctx[:xlimits]...)
        VTKView.yrange!(ctx[:plot],ctx[:flimits]...)
        VTKView.linewidth!(ctx[:plot],1)
    end
    if ctx[:clear]
        VTKView.clear!(ctx[:plot])
    end
    plot=ctx[:plot]
    VTKView.plotcolor!(plot,rgbtuple(ctx[:color])...)
    
    VTKView.title!(plot,ctx[:title])

    if ctx[:label]!=""
        VTKView.plotlegend!(plot,ctx[:label])
    end

    VTKView.addplot!(plot,collect(grid[Coordinates][1,:]),collect(func))
    reveal(ctx,TP)
end


