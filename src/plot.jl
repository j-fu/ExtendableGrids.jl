
"""
$(SIGNATURES)
Heuristically check if Plotter is VTKView
"""
isvtkview(Plotter)= (typeof(Plotter)==Module)&&isdefined(Plotter,:StaticFrame)

"""
$(SIGNATURES)
Heuristically check if Plotter is PyPlot
"""
ispyplot(Plotter)= (typeof(Plotter)==Module)&&isdefined(Plotter,:Gcf)

"""
$(SIGNATURES)
Heuristically check if  Plotter is Plots
"""
isplots(Plotter)= (typeof(Plotter)==Module) && isdefined(Plotter,:gr)


"""
$(SIGNATURES)
Heuristically check if Plotter is Makie/WGLMakie
"""
ismakie(Plotter)= (typeof(Plotter)==Module)&&isdefined(Plotter,:AbstractPlotting)

abstract type PyPlotType  end
abstract type MakieType   end
abstract type PlotsType   end
abstract type VTKViewType end

function plottertype(Plotter)
    if ismakie(Plotter)
        return MakieType
    elseif isplots(Plotter)
        return PlotsType
    elseif ispyplot(Plotter)
        return PyPlotType
    elseif isvtkview(Plotter)
        return VTKViewType
    end
    Nothing
end




"""
$(SIGNATURES)
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
    if ismakie(Plotter)
        return RGB(r,g,b)
    end
    if isplots(Plotter)
        return Plotter.RGB(r,g,b)
    end
end



const PlotterContext=Dict{Any,Any}
"""
$(SIGNATURES)
Initalize plotter context

Plotter: defaults to `nothing` and can be `PyPlot`, `Plots`, `VTKView`, `Makie`.

User controllable keys (depending on plot and plotter)
- :color:  color of plot on 1D grid
- :cmap:  color map for heatmap plot
- :label: label of plot
- :isolevels: number of isolevels
- :aspect: aspect ratio of plot
- :clear: if true (default) clear plot before plotting
- :show: if true (default) show plot
"""
PlotterContext(Plotter)=initialize_context(
    PlotterContext(:backend => plottertype(Plotter),
                   :Plotter => Plotter,
                   :colorlevels => 51,
                   :isolines => 11,
                   :colorbar => true,
                   :aspect => 1.0,
                   :show => true,
                   :axisgrid => true,
                   :clear => true,
                   :legend => true,
                   :legend_location => "upper right",
                   :colormap => "hot",
                   :label => "f",
                   :color => (0,0,0),
                   :elevation => false,
                   :elevation_factor => 1.0,
                   :resolution => (300,300)
                   ))

initialize_context(ctx::PlotterContext)=initialize_context(ctx,ctx[:backend])





plot!(ctx::PlotterContext,grid::ExtendableGrid)=plot!(ctx,ctx[:backend],Val{dim_space(grid)},grid)
plot!(ctx::PlotterContext,grid::ExtendableGrid, func)=plot!(ctx,ctx[:backend],Val{dim_space(grid)},grid,func)
plot!(ctx::PlotterContext,gf::GridFactory)=plot!(ctx,ctx[:backend],gf)


plot(grid::ExtendableGrid;Plotter=nothing)=plot!(PlotterContext(Plotter),grid)
plot(grid::ExtendableGrid,func ;Plotter=nothing)=plot!(PlotterContext(Plotter),grid,func)
plot(gf::GridFactory;Plotter=nothing)=plot!(PlotterContext(Plotter),gf)


initialize_context(ctx::PlotterContext,::Type{Nothing})=ctx
plot!(ctx, ::Type{Nothing}, ::Type{Val{1}},grid)=nothing
plot!(ctx, ::Type{Nothing}, ::Type{Val{2}},grid)=nothing
plot!(ctx, ::Type{Nothing}, ::Type{Val{3}},grid)=nothing

plot!(ctx, ::Type{Nothing}, ::Type{Val{1}},grid,func)=nothing
plot!(ctx, ::Type{Nothing}, ::Type{Val{2}},grid,func)=nothing
plot!(ctx, ::Type{Nothing}, ::Type{Val{3}},grid,func)=nothing


function bbox(grid)
    dim=dim_space(grid)
    coord=grid[Coordinates]
    bbmin=Vector(undef,dim)
    bbmax=Vector(undef,dim)
    bbmin.=1.0e30
    bbmax.=1.0e-30
    for i=1:size(coord,2)
        for idim=1:dim
            bbmin[idim]=min(bbmin[idim],coord[idim,i])
            bbmax[idim]=max(bbmax[idim],coord[idim,i])
        end
    end
    (bbmin, bbmax)
end
