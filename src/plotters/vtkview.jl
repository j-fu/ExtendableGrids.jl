function initialize_plot!(p::PlotContext,::Type{VTKViewType})
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

function plot!(ctx, ::Type{VTKViewType},grid)
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
    VTKView.display(frame)
end


plot!(ctx, T::Type{VTKViewType}, ::Type{Val{2}},grid)=plot!(ctx, T,grid)
plot!(ctx, T::Type{VTKViewType}, ::Type{Val{3}},grid)=plot!(ctx, T,grid)


function plot!(ctx, ::Type{VTKViewType},grid,func)
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
    VTKView.display(frame)
end



plot!(ctx, T::Type{VTKViewType}, ::Type{Val{2}},grid, func)=plot!(ctx, T,grid,func)
plot!(ctx, T::Type{VTKViewType}, ::Type{Val{3}},grid, func)=plot!(ctx, T,grid,func)
plot!(ctx, T::Type{VTKViewType}, ::Type{Val{1}},grid)=nothing

function plot!(ctx, T::Type{VTKViewType}, ::Type{Val{1}},grid, func)
    VTKView=ctx[:Plotter]
    frame=ctx[:frame]
    if !haskey(ctx,:plot)
        ctx[:plot]=VTKView.XYPlot()
        VTKView.addview!(frame,ctx[:plot],ctx[:iplot])
    end
    plot=ctx[:plot]
    VTKView.plotcolor!(plot,ctx[:color]...)
    VTKView.plotlegend!(plot,ctx[:label])
    VTKView.addplot!(plot,vec(grid[Coordinates]),func)
    display(frame)
end


