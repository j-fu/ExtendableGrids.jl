
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

"""
$(TYPEDEF)

Abstract type for dispatching on plotter
"""
abstract type PyPlotType  end

"""
$(TYPEDEF)

Abstract type for dispatching on plotter
"""
abstract type MakieType   end

"""
$(TYPEDEF)

Abstract type for dispatching on plotter
"""
abstract type PlotsType   end

"""
$(TYPEDEF)

Abstract type for dispatching on plotter
"""
abstract type VTKViewType end

"""
$(SIGNATURES)

Heuristically detect type of plotter, returns the corresponding 
abstract type fro plotting.
"""
function plottertype(Plotter::Union{Module,Nothing})
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


const SubPlotContext=Union{Dict{Symbol,Any},Nothing}
#
# Update plotter context from dict
#
function update_context!(ctx::SubPlotContext,kwargs)
    for (k,v) in kwargs
        ctx[Symbol(k)]=v
    end
    ctx
end



"""
$(TYPEDEF)

Context type for plots. Wrapper around dict.
"""
struct PlotterContext
    Plotter::Union{Module,Nothing}
    subplots::Array{SubPlotContext,2}
    context::SubPlotContext
    PlotterContext(Plotter::Union{Module,Nothing}, layout::Tuple, default::SubPlotContext)=new(Plotter, [copy(default) for I in CartesianIndices(layout)],copy(default))
end

layout(p::PlotterContext)=size(p.subplots)
subplot(p::PlotterContext, idx::Tuple)=p.subplots[idx...]
plottertype(p::PlotterContext)=plottertype(p.Plotter)
Base.copy(::Nothing)=nothing

#
# Default context information with help info.
#
default_plot_kwargs()=Dict{Symbol,Pair{Any,String}}(
    :colorlevels => Pair(51,"Number of color levels for contour plot"),
    :isolines => Pair(11,"Number of isolines in contour plot"),
    :colorbar => Pair(true,"Show colorbar in plots"),
    :aspect => Pair(1.0,"Aspect ration modification"),
    :show => Pair(true,"Show plots immediately"),
    :axisgrid => Pair(true,"Show background grid in plots"),
    :clear => Pair(true,"Clear plot before new plot."),
    :legend => Pair(true,"Add legend  to plot"),
    :legend_location => Pair("upper right","Location of legend"),
    :colormap => Pair("viridis","Contour plot colormap"),
    :label => Pair("f","Label of plot"),
    :layout => Pair((1,1),"Layout of plots"),
    :subplot => Pair((1,1),"Actual subplot"),
    :color => Pair((0,0,0),"Color of lines on plot"),
    :edges => Pair(true,"Plot grid edges when plotting grid"),
    :alpha => Pair(1.0,"Surface alpha value"),
    :interior => Pair(true,"Plot interior of grid"),
    :outline => Pair(true,"Plot outline of domain"),
    :xplane => Pair(prevfloat(Inf),"xplane for 3D visualization"),
    :yplane => Pair(prevfloat(Inf),"yplane for 3D visualization"),
    :zplane => Pair(prevfloat(Inf),"zplane for 3D visualization"),
    :flevel => Pair(prevfloat(Inf),"isolevel for 3D visualization"),
    :azim => Pair(-60,"azimuth angle for 3D visualization (in degrees)"),
    :elev => Pair(30,"Elevation angle for 3D visualization (in degrees)"),
    :title => Pair(" ","Plot title"),
    :elevation => Pair(0.0,"Height factor for elevation of 2D plot"),
    :resolution => Pair((500,500),"Plot xy resolution"),
    :framepos => Pair(1,"Subplot position in frame (VTKView)"),
    :fignumber => Pair(1,"Figure number (PyPlot)")
)

#
# Print default dict for interpolation into docstrings
#
function myprint(dict::Dict{Symbol,Pair{Any,String}})
    lines_out=IOBuffer()
    for (k,v) in dict
        println(lines_out,"  - $(k): $(v[2]). Default: $(v[1])\n")
    end
    String(take!(lines_out))
end

"""
$(SIGNATURES)

Initalize plotter context

Plotter: defaults to `nothing` and can be `PyPlot`, `Plots`, `VTKView`, `Makie`.

Keyword arguments:

$(myprint(default_plot_kwargs()))
"""
function PlotterContext(;Plotter::Union{Module,Nothing}=nothing, kwargs...)
    default_ctx=Dict{Symbol,Any}( k => v[1] for (k,v) in default_plot_kwargs())
    update_context!(default_ctx,kwargs)
    layout=default_ctx[:layout]
    if isnothing(Plotter)
        default_ctx=nothing
    end
    p=PlotterContext(Plotter,layout,default_ctx)
    if !isnothing(Plotter)
        p.context[:Plotter]=Plotter
        initialize_plot!(p,plottertype(Plotter))
        for I in CartesianIndices(layout)
            ctx=p.subplots[I]
            i=Tuple(I)
            ctx[:subplot]=i
            ctx[:iplot]=layout[2]*(i[1]-1)+i[2]
            ctx[:Plotter]=Plotter
        end
    end
    p
end

Base.getindex(p::PlotterContext,i,j)=p.subplots[i,j]

"""
$(SIGNATURES)

Plot grid.

Keyword arguments:

$(myprint(default_plot_kwargs()))
"""
plot!(ctx::SubPlotContext,grid::ExtendableGrid; kwargs...)=plot!(update_context!(ctx,kwargs),plottertype(ctx[:Plotter]),Val{dim_space(grid)},grid)
plot!(p::PlotterContext,grid::ExtendableGrid; kwargs...)=plot!(p[1,1],grid,kwargs...)

"""
$(SIGNATURES)

Plot grid without predefined context.

Keyword arguments:

$(myprint(default_plot_kwargs()))
"""
plot(grid::ExtendableGrid; Plotter=nothing, kwargs...)=plot!(PlotterContext(Plotter=Plotter; kwargs...),grid)



"""
$(SIGNATURES)

Plot vector on grid as P1 FEM function.

Keyword arguments

$(myprint(default_plot_kwargs()))
"""
function plot!(ctx::SubPlotContext,grid::ExtendableGrid,func; kwargs...)
    plot!(update_context!(ctx,kwargs),plottertype(ctx[:Plotter]),Val{dim_space(grid)},grid,func)
end

plot!(p::PlotterContext, grid::ExtendableGrid, func; kwargs...)=plot!(p[1,1],grid,func,kwargs...)


"""
$(SIGNATURES)

Plot vector on grid without predefined context as P1 FEM function.

Keyword arguments:

$(myprint(default_plot_kwargs()))
"""
plot(grid::ExtendableGrid,func ;Plotter=nothing,kwargs...)=plot!(PlotterContext(Plotter=Plotter;kwargs...),grid,func)




#
# Dummy functions to allow Plotter=nothing
#
initialize_subplot(ctx::SubPlotContext,::Type{Nothing})=ctx
update_context!(::Nothing,kwargs)=nothing

plot!(ctx::Nothing,grid::ExtendableGrid;kwargs...)=nothing
plot!(ctx::Nothing,grid::ExtendableGrid,func;kwargs...)=nothing

plot!(ctx, ::Type{Nothing}, ::Type{Val{1}},grid)=nothing
plot!(ctx, ::Type{Nothing}, ::Type{Val{2}},grid)=nothing
plot!(ctx, ::Type{Nothing}, ::Type{Val{3}},grid)=nothing

plot!(ctx, ::Type{Nothing}, ::Type{Val{1}},grid,func)=nothing
plot!(ctx, ::Type{Nothing}, ::Type{Val{2}},grid,func)=nothing
plot!(ctx, ::Type{Nothing}, ::Type{Val{3}},grid,func)=nothing

