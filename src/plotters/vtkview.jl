function initialize_context!(ctx::PlotterContext,::Type{VTKViewType})
    ctx
end

function prepare_frame!(ctx)
    VTKView=ctx[:Plotter]
    if !haskey(ctx,:frame)
        ctx[:frame]=VTKView.StaticFrame()
    end
    VTKView.size!(ctx[:frame],ctx[:resolution]...)
    ctx
end

function plot!(ctx, ::Type{VTKViewType},grid)
    VTKView=ctx[:Plotter]
    prepare_frame!(ctx)
    frame=ctx[:frame]
    VTKView.clear!(frame)
    dataset=VTKView.DataSet()
    VTKView.simplexgrid!(dataset,grid[Coordinates],grid[CellNodes])
    VTKView.boundarygrid!(dataset,grid[BFaceNodes])
    VTKView.boundarymarker!(dataset,grid[BFaceRegions])
    VTKView.cellmarker!(dataset,grid[CellRegions])
    gridview=VTKView.GridView()
    VTKView.data!(gridview,dataset)
    VTKView.addview!(frame,gridview)
    VTKView.display(frame)
end


plot!(ctx, T::Type{VTKViewType}, ::Type{Val{2}},grid)=plot!(ctx, T,grid)
plot!(ctx, T::Type{VTKViewType}, ::Type{Val{3}},grid)=plot!(ctx, T,grid)


function plot!(ctx, ::Type{VTKViewType},grid,func)
    VTKView=ctx[:Plotter]
    prepare_frame!(ctx)
    frame=ctx[:frame]
    if !haskey(ctx,:dataset) || !seemingly_equal(grid,ctx[:grid])
        VTKView.clear!(frame)
        ctx[:grid]=grid
        ctx[:dataset]=VTKView.DataSet()
        VTKView.simplexgrid!(ctx[:dataset],grid[Coordinates],grid[CellNodes])
        scalarview=VTKView.ScalarView()
        VTKView.data!(scalarview,ctx[:dataset],ctx[:label])
        VTKView.addview!(frame,scalarview)
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
        VTKView.clear!(frame)
        ctx[:plot]=VTKView.XYPlot()
        VTKView.addview!(frame,ctx[:plot])
    end
    plot=ctx[:plot]
    VTKView.plotcolor!(plot,ctx[:color]...)
    VTKView.plotlegend!(plot,ctx[:label])
    VTKView.addplot!(plot,vec(grid[Coordinates]),func)
    display(frame)
end


