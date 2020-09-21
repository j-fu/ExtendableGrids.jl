function initialize_context!(ctx::PlotterContext,::Type{VTKViewType})
    ctx
end

function prepare_frame!(ctx)
    VTKView=ctx[:Plotter]
    if !haskey(ctx,:frame)
        ctx[:frame]=VTKView.StaticFrame()
        VTKView.layout!(ctx[:frame],ctx[:layout]...)
        VTKView.size!(ctx[:frame],ctx[:resolution]...)
    end
    ctx
end

function plot!(ctx, ::Type{VTKViewType},grid)
    VTKView=ctx[:Plotter]
    prepare_frame!(ctx)
    frame=ctx[:frame]
    if !haskey(ctx,:dataset)
        if ctx[:clear]
            VTKView.clear!(frame)
        end
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
        VTKView.addview!(frame,ctx[:gridview],ctx[:framepos])
    end
    VTKView.display(frame)
end


plot!(ctx, T::Type{VTKViewType}, ::Type{Val{2}},grid)=plot!(ctx, T,grid)
plot!(ctx, T::Type{VTKViewType}, ::Type{Val{3}},grid)=plot!(ctx, T,grid)


function plot!(ctx, ::Type{VTKViewType},grid,func)
    VTKView=ctx[:Plotter]
    prepare_frame!(ctx)
    frame=ctx[:frame]
    if !haskey(ctx,:dataset)
        if ctx[:clear]
            VTKView.clear!(frame)
        end
        ctx[:dataset]=VTKView.DataSet()
    end
    if !haskey(ctx,:grid)
        ctx[:grid]=grid
        VTKView.simplexgrid!(ctx[:dataset],grid[Coordinates],grid[CellNodes])
    end
    if !haskey(ctx,:scalarview)
        ctx[:scalarview]=VTKView.ScalarView()
        VTKView.addview!(frame,ctx[:scalarview],ctx[:framepos])
        VTKView.data!(ctx[:scalarview],ctx[:dataset],ctx[:label])
    end
    if !seemingly_equal(grid,ctx[:grid])
        VTKView.simplexgrid!(ctx[:dataset],grid[Coordinates],grid[CellNodes])
    end
    VTKView.pointscalar!(ctx[:dataset],func,ctx[:label])
    VTKView.display(frame)
end



plot!(ctx, T::Type{VTKViewType}, ::Type{Val{2}},grid, func)=plot!(ctx, T,grid,func)
plot!(ctx, T::Type{VTKViewType}, ::Type{Val{3}},grid, func)=plot!(ctx, T,grid,func)


function plot!(ctx, T::Type{VTKViewType}, ::Type{Val{1}},grid, func)
    VTKView=ctx[:Plotter]
    prepare_frame!(ctx)
    frame=ctx[:frame]
    if !haskey(ctx,:plot)
        if ctx[:clear]
            VTKView.clear!(frame)
        end
        ctx[:plot]=VTKView.XYPlot()
        VTKView.addview!(frame,ctx[:plot],ctx[:framepos])
    end
    plot=ctx[:plot]
    VTKView.plotcolor!(plot,ctx[:color]...)
    VTKView.plotlegend!(plot,ctx[:label])
    VTKView.addplot!(plot,vec(grid[Coordinates]),func)
    display(frame)
end


