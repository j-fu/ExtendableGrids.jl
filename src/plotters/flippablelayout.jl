module FlippableLayout

# Thanks to Julius Krummbiegel for providing
# a basic implementation to focus switching


# Currently this module sits within the plotting section of ExtendableGrids.jl
# which avoids creating dependencies on plotting backends.
# So we provide a way to emulate "import Makie" by allowing
# to set it as a global variable.

Makie=nothing

"""
   FLayout

Struct describing flippable layout data

"""
mutable struct FLayout
    visible #::GridLayout
    offscreen #::GridLayout
    blocked::Bool
    layoutables::Dict{Tuple{Int64,Int64},Any} # Any -> Layoutable
    condition::Condition
    FLayout(visible;blocked=false)=new(visible,
                                       Makie.GridLayout(bbox = Makie.BBox(-500, -400, -500, -400)),
                                       blocked,
                                       Dict{Tuple{Int64,Int64},Any}(),
                                       Condition())
end


function Base.setindex!(flayout::FLayout,layoutable,i,j)
    if isnothing(layoutable)
        flayout.offscreen[1, 1] = flayout.layoutables[(i,j)]
#        delete!(flayout.visible,flayout.layoutables[(i,j)])#may be this does not work anymore
        delete!(flayout.layoutables,(i,j)) 
    elseif !isa(layoutable,Makie.MakieLayout.Layoutable)
        error("can only set layoutables")
    else
        flayout.layoutables[(i,j)]=layoutable
        _showall(flayout)
        yield()
    end
end


Base.getindex(flayout::FLayout,i,j)=flayout.layoutables[(i,j)]


function _showall(flayout::FLayout)
    for (pos, layoutable) in flayout.layoutables
        flayout.visible[pos...]=layoutable
    end
    Makie.trim!(flayout.visible)
end

#
# Check if mouse position is  within pixel area of scene
#
function _inscene(scene,pos)
    area=scene.px_area[]
    pos[1]>area.origin[1] &&
        pos[1] < area.origin[1]+area.widths[1] &&
        pos[2]>area.origin[2] &&
        pos[2] < area.origin[2]+area.widths[2]
end


"""
   flayoutscene(;blocked=false, kwargs...)

Layoutscene with interactive layout and blocking functionality

The `,` key switches between focused view showing only one subscene
and "gallery view" showing all layoutables at once.

The space key toggles blocking of the execution of the main therad
when `yield` is replaced by `yieldwait`. Initial blocking state is 
set by the `blocked` kwarg.

The `kwargs...` are the same as of `AbstractPlotting.layoutscene`.


The idea is that this can work in some cases as a drop-in replacement
of `layoutscene`.     
"""

function flayoutscene(;blocked=false,
                      focuskey=Makie.Keyboard.comma,
                      blockingkey=Makie.Keyboard.space,
                      kwargs...)
    
    (parent, layout) = Makie.layoutscene(;kwargs...)
    
    flayout=FLayout(layout,blocked=blocked)

    gallery_view=Makie.Node(true)

    # Watch mouse position
    mouseposition=Makie.Node((0,0))

    Makie.on(m->mouseposition[]=m, parent.events.mouseposition)

    # Switch focus to subscene  at pos
    function _focus(focus)
        for (key,layoutable) in flayout.layoutables
            if key==focus
                flayout.visible[1, 1] = layoutable
            else
                flayout.offscreen[1, 1] = layoutable
            end
        end
        Makie.trim!(flayout.visible)
    end
    

    # Figure out to which subscene the mouse position
    # corresponds
    function _subscene(mouseposition)
        for (key,layoutable) in flayout.layoutables
            if _inscene(layoutable.scene,mouseposition)
                return key
            end
        end
        return nothing
    end
    
    function _toggle_block(flayout::FLayout)
        if flayout.blocked
            flayout.blocked=false
            notify(flayout.condition)
        else
            flayout.blocked=true
        end
    end
    
    # Handle global key events for `,` (focus/gallery view)
    # and space (toggle blocking)
    Makie.on(parent.events.keyboardbuttons) do buttons
        if Makie.ispressed(parent, focuskey)
            if gallery_view[]
                pos=_subscene(mouseposition[])
                if !isnothing(pos)
                    _focus(pos)
                end
                gallery_view[]=false
            else
                _showall(flayout)
                gallery_view[]=true
            end
        end
        if Makie.ispressed(parent, blockingkey)
            _toggle_block(flayout)
        end
    end
    parent,flayout
end

"""
     yieldwait(fliplayoutscene)

Yield and wait in case of scene being blocked via space key toggle
"""
function yieldwait(flayout::FLayout)
    yield()
    if flayout.blocked
        wait(flayout.condition)
    end
end

"""
    setmakie!(MyMakie)

Set the Makie module.
"""
setmakie!(MyMakie) = global Makie=MyMakie


export flayoutscene
export yieldwait
export setmakie!

end
