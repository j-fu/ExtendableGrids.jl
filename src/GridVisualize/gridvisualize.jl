module GridVisualize

using ExtendableGrids
using DocStringExtensions
using ElasticArrays
using StaticArrays

using Colors
using ColorSchemes
using GeometryBasics
using LinearAlgebra


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
$(SIGNATURES)
Heuristically check if Plotter is MeshCat
"""
ismeshcat(Plotter)= (typeof(Plotter)==Module)&&isdefined(Plotter,:Visualizer)

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
$(TYPEDEF)

Abstract type for dispatching on plotter
"""
abstract type MeshCatType end

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
    elseif ismeshcat(Plotter)
        return MeshCatType
    end
    Nothing
end


"""
$(TYPEDEF)

A SubVis is just a dictionary which contains plotting information,
including type of the plotter and its position in the plot.
"""
const SubVis=Union{Dict{Symbol,Any},Nothing}

#
# Update subplot context from dict
#
function _update_context!(ctx::SubVis,kwargs)
    for (k,v) in kwargs
        ctx[Symbol(k)]=v
    end
    ctx
end



"""
$(TYPEDEF)

Context type for plots.
"""
struct GridVisualizer
    Plotter::Union{Module,Nothing}
    subplots::Array{SubVis,2}
    context::SubVis
    GridVisualizer(Plotter::Union{Module,Nothing}, layout::Tuple, default::SubVis)=new(Plotter,
                                                                                            [copy(default) for I in CartesianIndices(layout)],
                                                                                            copy(default))
end


"""
$(SIGNATURES)

Return the layout of a GridVisualizer
"""
Base.size(p::GridVisualizer)=size(p.subplots)

"""
$(SIGNATURES)

Return a SubVis
"""
Base.getindex(p::GridVisualizer,i,j)=p.subplots[i,j]


"""
$(SIGNATURES)

Return the type of a plotter.
"""
plottertype(p::GridVisualizer)=plottertype(p.Plotter)

#
# Default context information with help info.
#
default_plot_kwargs()=Dict{Symbol,Pair{Any,String}}(
    :colorlevels => Pair(51,"Number of color levels for contour plot"),
    :isolines => Pair(11,"Number of isolines in contour plot"),
    :colorbar => Pair(true,"Show colorbar in plots"),
    :aspect => Pair(1.0,"Aspect ratio modification"),
    :show => Pair(false,"Show plot immediately"),
    :reveal => Pair(false,"Show plot immediately (same as :show)"),
    :fast => Pair(true,"Fast updates (with Makie, restricted functionality)"),
    :cellwise => Pair(false,"Cellwise 1D plot, can be slow)"),
    :axisgrid => Pair(true,"Show background grid in plots"),
    :clear => Pair(true,"Clear plot before new plot."),
    :legend => Pair(true,"Add legend  to plot"),
    :legend_location => Pair("upper right","Location of legend"),
    :colormap => Pair(:viridis,"Contour plot colormap (any from [ColorSchemes.jl](https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#Pre-defined-schemes))"),
    :label => Pair("","Label of plot"),
    :xlimits => Pair((1,-1),"x limits"),
    :ylimits => Pair((1,-1),"y limits"),
    :zlimits => Pair((1,-1),"z limits"),
    :flimits => Pair((1,-1),"function limits"),
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
    :azim => Pair(-60,"Azimuth angle for 3D visualization (in degrees)"),
    :elev => Pair(30,"Elevation angle for 3D visualization (in degrees)"),
    :title => Pair("","Plot title"),
    :elevation => Pair(0.0,"Height factor for elevation of 2D plot"),
    :resolution => Pair((500,500),"Plot xy resolution"),
    :framepos => Pair(1,"Subplot position in frame (VTKView)"),
    :fignumber => Pair(1,"Figure number (PyPlot)")
)

#
# Print default dict for interpolation into docstrings
#
function _myprint(dict::Dict{Symbol,Pair{Any,String}})
    lines_out=IOBuffer()
    for (k,v) in dict
        println(lines_out,"  - $(k): $(v[2]). Default: $(v[1])\n")
    end
    String(take!(lines_out))
end

"""
$(SIGNATURES)

Create a plot context.

Plotter: defaults to `nothing` and can be `PyPlot`, `Plots`, `VTKView`, `Makie`.
This pattern to pass the backend as a module to a plot function allows to circumvent
to create heavy default package dependencies.


Depending on the `layout` keyword argument, a 2D grid of subplots is created.
Further `visualize!` commands then plot into one of these subplots:

```julia
p=GridVisualizer(Plotter=PyPlot, layout=(2,2)
visualize!(p[1,2], ...)
````

A `plot`  command just implicitely creates a plot context:
```julia
visualize(..., Plotter=PyPlot) 
```
is equivalent to
```julia
p=GridVisualizer(Plotter=PyPlot, layout=(1,1)
visualize!(p,...) 
```

Please not that the return values of all plot commands are specific to the Plotter.

Depending on the backend, interactive mode switch between "gallery view" showing all plots at
onece and "focused view" showing only one plot is possible.


Keyword arguments:

$(_myprint(default_plot_kwargs()))
"""
function GridVisualizer(;Plotter::Union{Module,Nothing}=nothing, kwargs...)
    default_ctx=Dict{Symbol,Any}( k => v[1] for (k,v) in default_plot_kwargs())
    _update_context!(default_ctx,kwargs)
    layout=default_ctx[:layout]
    if isnothing(Plotter)
        default_ctx=nothing
    end
    p=GridVisualizer(Plotter,layout,default_ctx)
    if !isnothing(Plotter)
        p.context[:Plotter]=Plotter
        for I in CartesianIndices(layout)
            ctx=p.subplots[I]
            i=Tuple(I)
            ctx[:subplot]=i
            ctx[:iplot]=layout[2]*(i[1]-1)+i[2]
            ctx[:Plotter]=Plotter
            ctx[:GridVisualizer]=p
        end
        initialize!(p,plottertype(Plotter))
    end
    p
end


"""
$(SIGNATURES)

Plot grid.

Keyword arguments:

$(_myprint(default_plot_kwargs()))
"""
function visualize!(ctx::SubVis,grid::ExtendableGrid; kwargs...)
    _update_context!(ctx,kwargs)
    visualize!(ctx,plottertype(ctx[:Plotter]),Val{dim_space(grid)},grid)
end

"""
$(SIGNATURES)

Plot grid.

Keyword arguments:

$(_myprint(default_plot_kwargs()))
"""
visualize!(p::GridVisualizer,grid::ExtendableGrid; kwargs...)=visualize!(p[1,1],grid,kwargs...)


"""
$(SIGNATURES)

Plot grid without predefined context.

Keyword arguments:

$(_myprint(default_plot_kwargs()))
"""
visualize(grid::ExtendableGrid; Plotter=nothing, kwargs...)=visualize!(GridVisualizer(Plotter=Plotter; show=true, kwargs...),grid)


"""
$(SIGNATURES)

Plot vector on grid as P1 FEM function.

Keyword arguments

$(_myprint(default_plot_kwargs()))
"""
function visualize!(ctx::SubVis,grid::ExtendableGrid,func; kwargs...)
    _update_context!(ctx,Dict(:clear=>true,:show=>false,:reveal=>false))
    _update_context!(ctx,kwargs)
    visualize!(ctx,plottertype(ctx[:Plotter]),Val{dim_space(grid)},grid,func)
end


"""
$(SIGNATURES)

Plot vector on grid without predefined context as P1 FEM function.

Keyword arguments:

$(_myprint(default_plot_kwargs()))
"""
visualize(grid::ExtendableGrid,func ;Plotter=nothing,kwargs...)=visualize!(GridVisualizer(Plotter=Plotter;kwargs...),grid,func,show=true)

visualize!(p::GridVisualizer,grid::ExtendableGrid, func; kwargs...)=visualize!(p[1,1],grid,func; kwargs...)
visualize!(p::GridVisualizer,grid::ExtendableGrid, kwargs...)=visualize!(p[1,1],grid; kwargs...)


visualize(X::Vector,func ;kwargs...)=visualize(simplexgrid(X),func;kwargs...)
visualize!(ctx::SubVis,X::Vector,func; kwargs...)=visualize!(ctx,simplexgrid(X),func;kwargs...)
visualize!(ctx::GridVisualizer,X::Vector,func; kwargs...)=visualize!(ctx,simplexgrid(X),func;kwargs...)

"""
$(SIGNATURES)

Finish and show plot. Same as setting `:reveal=true` or `:show=true` in last visualize statment
for a context.
"""
reveal(p::GridVisualizer)=reveal(p, plottertype(p.Plotter))

"""
$(SIGNATURES)

Save figure to disk
"""
save(fname,p::GridVisualizer)=save(fname,p, plottertype(p.Plotter))


#
# Dummy methods to allow Plotter=nothing
#
_update_context!(::Nothing,kwargs)=nothing
Base.copy(::Nothing)=nothing
visualize!(ctx::Nothing,grid::ExtendableGrid;kwargs...)=nothing
visualize!(ctx::Nothing,grid::ExtendableGrid,func;kwargs...)=nothing

visualize!(ctx, ::Type{Nothing}, ::Type{Val{1}},grid)=nothing
visualize!(ctx, ::Type{Nothing}, ::Type{Val{2}},grid)=nothing
visualize!(ctx, ::Type{Nothing}, ::Type{Val{3}},grid)=nothing

visualize!(ctx, ::Type{Nothing}, ::Type{Val{1}},grid,func)=nothing
visualize!(ctx, ::Type{Nothing}, ::Type{Val{2}},grid,func)=nothing
visualize!(ctx, ::Type{Nothing}, ::Type{Val{3}},grid,func)=nothing


displayable(ctx,Any)=nothing
reveal(p,::Type{Nothing})=nothing





include("common.jl")
include("pyplot.jl")
include("makie.jl")
include("vtkview.jl")
include("meshcat.jl")
include("plots.jl")


export visualize,visualize!,save,reveal
export isplots,isvtkview,ispyplot,ismakie
export GridVisualizer, SubVis
export plottertype
export PyPlotType,MakieType,PlotsType,VTKViewType 

end

