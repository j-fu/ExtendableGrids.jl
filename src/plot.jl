
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
$(TYPEDEF)

Context type for plots. Wrapper around dict.
"""
struct PlotterContext<: AbstractDict{Any,Any}
    dict::Dict{Any,Any}
    PlotterContext(default_values::Dict{Any,Any})=new(default_values)
end

Base.setindex!(ctx::PlotterContext,v,k)=setindex!(ctx.dict,v,k)
Base.getindex(ctx::PlotterContext,k)=getindex(ctx.dict,k)
Base.keys(ctx::PlotterContext)=keys(ctx.dict)
Base.haskey(ctx::PlotterContext)=haskey(ctx.dict)
Base.iterate(ctx::PlotterContext) = iterate(ctx.dict)
Base.iterate(ctx::PlotterContext,state) = iterate(ctx.dict,state)
Base.length(ctx::PlotterContext)=length(ctx.dict)

#
# Default context information with help info.
#
default_plot_kwargs()=Dict{Any,Pair{Any,String}}(
    :colorlevels => Pair(51,"Number of color levels for contour plot"),
    :isolines => Pair(11,"Number of isolines in contour plot"),
    :colorbar => Pair(true,"Show colorbar in plots"),
    :aspect => Pair(1.0,"Aspect ration modification"),
    :show => Pair(true,"Show plots immediately"),
    :axisgrid => Pair(true,"Show background grid in plots"),
    :clear => Pair(true,"Clear plot before new plot."),
    :legend => Pair(true,"Add legend  to plot"),
    :legend_location => Pair("upper right","Location of legend"),
    :colormap => Pair("hot","Contour plot colormap"),
    :label => Pair("f","Label of plot"),
    :layout => Pair((1,1),"Layout of multiple plots (experimental)"),
    :subplot => Pair(111,"Subplot number for multiple plots"),
    :color => Pair((0,0,0),"Color of lines on plot"),
    :edges => Pair(true,"Plot grid edges when plotting grid"),
    :alpha => Pair(1.0,"Surface alpha value"),
    :interior => Pair(true,"plot interior of grid"),
    :xplane => Pair(1.0,"xplane for 3D visualization"),
    :yplane => Pair(1.0,"yplane for 3D visualization"),
    :zplane => Pair(1.0,"zplane for 3D visualization"),
    :flevel => Pair(1.0,"isolevel for 3D visualization"),
    :azim => Pair(-60,"azimuth angle for 3D visualization (in degrees)"),
    :elev => Pair(30,"elevation angle for 3D visualization (in degrees)"),
    :elevation => Pair(0.0,"Height factor for elevation of 2D plot"),
    :resolution => Pair((500,500),"Plot xy resolution"),
    :framepos => Pair(1,"Subplot position in frame (VTKView)"),
    :fignumber => Pair(1,"Figure number (PyPlot)")
)

#
# Print default dict for interpolation into docstrings
#
function myprint(dict::Dict{Any,Pair{Any,String}})
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
function PlotterContext(Plotter::Module; kwargs...)
    ctx=PlotterContext(Dict{Any,Any}(:backend => plottertype(Plotter),
                                     :Plotter => Plotter,
                                     ))
    update_context!(ctx,Dict{Any,Any}( k => v[1] for (k,v) in default_plot_kwargs()))
    initialize_context!(ctx,ctx[:backend])
    update_context!(ctx,kwargs)
end


#
# Update plotter context from dict
#
function update_context!(ctx::PlotterContext,kwargs)
    for (k,v) in kwargs
        ctx[Symbol(k)]=v
    end
    ctx
end


"""
$(SIGNATURES)

Plot grid.

Keyword arguments:

$(myprint(default_plot_kwargs()))
"""
plot!(ctx::PlotterContext,grid::ExtendableGrid;kwargs...)=plot!(update_context!(ctx,kwargs),ctx[:backend],Val{dim_space(grid)},grid)

"""
$(SIGNATURES)

Plot vector on grid as P1 FEM function.

Keyword arguments:

$(myprint(default_plot_kwargs()))
"""
plot!(ctx::PlotterContext,grid::ExtendableGrid,func;kwargs...)=plot!(update_context!(ctx,kwargs),ctx[:backend],Val{dim_space(grid)},grid,func)


"""
$(SIGNATURES)

Plot grid without predefined context.

Keyword arguments:

$(myprint(default_plot_kwargs()))
"""
plot(grid::ExtendableGrid;Plotter=nothing, kwargs...)=plot!(update_context!(PlotterContext(Plotter),kwargs),grid)


"""
$(SIGNATURES)

Plot vector on grid without predefined context as P1 FEM function.

Keyword arguments:

$(myprint(default_plot_kwargs()))
"""
plot(grid::ExtendableGrid,func ;Plotter=nothing,kwargs...)=plot!(update_context!(PlotterContext(Plotter),kwargs),grid,func)




#
# Dummy functions to allow Plotter=nothing
#
PlotterContext(::Nothing; kwargs...)=nothing
initialize_context(ctx::PlotterContext,::Type{Nothing})=ctx
update_context!(::Nothing,kwargs)=nothing

plot!(ctx::Nothing,grid::ExtendableGrid;kwargs...)=nothing
plot!(ctx::Nothing,grid::ExtendableGrid,func;kwargs...)=nothing

plot!(ctx, ::Type{Nothing}, ::Type{Val{1}},grid)=nothing
plot!(ctx, ::Type{Nothing}, ::Type{Val{2}},grid)=nothing
plot!(ctx, ::Type{Nothing}, ::Type{Val{3}},grid)=nothing

plot!(ctx, ::Type{Nothing}, ::Type{Val{1}},grid,func)=nothing
plot!(ctx, ::Type{Nothing}, ::Type{Val{2}},grid,func)=nothing
plot!(ctx, ::Type{Nothing}, ::Type{Val{3}},grid,func)=nothing

