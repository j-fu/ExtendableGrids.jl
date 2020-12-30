module MultiScene

# Thanks to Julius Krummbiegel for providing
# a basic implementation to focus switching


# Currently this module sits within the plotting section of ExtendableGrids.jl
# which avoids creating dependencies on plotting backends.
# So we provide a way to emulate "import Makie" by allowing
# to set it as a global variable.

Makie=nothing



#
# Check if position is  within pixel area of scene
#
function _inscene(scene,pos)
    area=scene.px_area[]
    pos[1]>area.origin[1] &&
        pos[1] < area.origin[1]+area.widths[1] &&
        pos[2]>area.origin[2] &&
        pos[2] < area.origin[2]+area.widths[2]
end


# Calculate centered text position
_textposition(scene)=Makie.lift(a->Makie.Vec2f0(Makie.widths(a)[1] ./ 2 , 0), Makie.pixelarea(scene))

# Create header text
_header(scene,txt)=Makie.text(txt,
                              textsize = 20,
                              raw=true,
                              position=_textposition(scene),
                              camera=Makie.campixel!,
                              align = (:center, :bottom))

# Create footer text
_footer(scene,status)=Makie.text(Makie.lift(a->a,status),
                                 textsize = 15,raw=true,
                                 position=_textposition(scene),
                                 camera=Makie.campixel!,
                                 align = (:center, :bottom))



# Handle blocking/unblocking via space key
mutable struct Blocker
    condition::Condition
    blocked::Bool
    Blocker(;blocked=false)=new(Condition(),blocked)
end


# block/unblock, and notify about unblocking
function _toggle(blocker::Blocker)
    if blocker.blocked
        blocker.blocked=false
        notify(blocker.condition)
    else
        blocker.blocked=true
    end
end

# yield and wait for block 
function yieldwait(blocker::Blocker)
    yield()
    if blocker.blocked
        wait(blocker.condition)
    end
end

#
# As we have only one makie window, we also can afford
# to have one  blocker.
#
globalblocker=Blocker()



"""
     multiscene(;layout=(1,1), blocked=false, resolution=(500,500))


Create a scene with given layout grid and resolution  of the 
window.  Returns an array of subscenes according to the layout parameter.

The `,` key switches between focused view showing only one subscene
and "gallery view" showing all subscenes at once.

The space key toggles blocking of the execution of the main therad
when `yield` is replaced by `yieldwait`.

"""
function multiscene(;layout=(1,1), blocked=false, resolution=(500,500))

    (scene, makielayout) = Makie.layoutscene(resolution = resolution)

    offscreen_gl = Makie.GridLayout(bbox = Makie.BBox(-500, -400, -500, -400))

    axs = makielayout[] = [Makie.LScene(scene) for _ in CartesianIndices(layout)]

    gallery_view=Makie.Node(true)

    globalblocker.blocked=blocked

    # Watch mouse position
    mouseposition=Makie.Node((0,0))
    Makie.on(m->mouseposition[]=m, scene.events.mouseposition)

    # Switch focus to i-th subscene 
    function _focus(i)
        for (j, ax) in enumerate(axs)
            if j != i
                offscreen_gl[1, 1] = ax
            else
                makielayout[1, 1] = ax
            end
        end
        Makie.trim!(makielayout)
        gallery_view[]=false
    end

    # Switch back to gallery view
    function _show_all()
        for idx in CartesianIndices(axs)
            makielayout[Tuple(idx)...]=  axs[Tuple(idx)...]
        end
        gallery_view[]=true
    end

    # Figure out to which subscene the mouse position
    # corresponds
    function _subscene(mouseposition)
        for i=1:length(scene.children)
            if _inscene(scene.children[i],mouseposition)
                return i
            end
        end
        return 0
    end
    
    # Handle global key events for `,` (focus/gallery view)
    # and space (toggle blocking)
    Makie.on(scene.events.keyboardbuttons) do buttons
        if Makie.ispressed(scene, Makie.Keyboard.comma)
            gallery_view[] ? _focus(_subscene(mouseposition[])) : _show_all() 
        end
        if Makie.ispressed(scene, Makie.Keyboard.space)
            _toggle(globalblocker)
        end
    end
    scene,axs
end


"""
     add_scene(ax, scene, title=" ",  staus=nothing)

Combine scene with title header and status footer and add it to ax.
The optional staus footer should be a Node containing text.
"""
function add_scene!(ax,scene;title=" ",status=nothing)
    if isnothing(status)
        Makie.hbox(scene,_header(scene,title),parent=ax.scene)
    else
        Makie.hbox(_footer(scene,status),scene,_header(scene,title),parent=ax.scene)
    end
end


"""
     scene_interaction(update_scene,scene,switchkeys::Vector{Symbol}=[:nothing])   

Control multiple scene elements via keyboard up/down keys. 
Each switchkey is assumed to correspond to one of these elements.
Pressing a switch key transfers control to its associated element.

Control of values of the current associated element is performed
by triggering change values via up/down (± 1)  resp. page_up/page_down (±10) keys

The update_scene callbac gets passed the change value and the symbol.
"""
function scene_interaction(update_scene,scene,switchkeys::Vector{Symbol}=[:nothing])

    # Initial active switch key is the first in the vector passed
    activeswitch=Makie.Node(switchkeys[1])

    # Handle mouse position within scene
    mouseposition=Makie.Node((0,0))
    Makie.on(m->mouseposition[]=m, scene.events.mouseposition)

    # Set keyboard event callback
    Makie.on(scene.events.keyboardbuttons) do buttons
        if _inscene(scene,mouseposition[])
            # On pressing a switch key, pass control
            for i=1:length(switchkeys)
                if switchkeys[i]!=:nothing && Makie.ispressed(scene,getproperty(Makie.Keyboard,switchkeys[i]))
                    activeswitch[]=switchkeys[i]
                    update_scene(0,switchkeys[i])
                    return
                end
            end
            
            # Handle change values via up/down control
            if Makie.ispressed(scene, Makie.Keyboard.up)
                update_scene(1,activeswitch[])
            elseif Makie.ispressed(scene, Makie.Keyboard.down)
                update_scene(-1,activeswitch[])
            elseif Makie.ispressed(scene, Makie.Keyboard.page_up)
                update_scene(10,activeswitch[])
            elseif Makie.ispressed(scene, Makie.Keyboard.page_down)
                update_scene(-10,activeswitch[])
            end
        end
    end
end



"""
     yieldwait()

Yield and wait in case of scene being blocked via space key toggle
"""
yieldwait()=yieldwait(globalblocker)

"""
      subscenecontext(ax; ctx=Dict{Symbol,Any}())

    Set/Get subscene context as attribute to ax
"""
function subscenecontext(ax; ctx=Dict{Symbol,Any}())
    if !haskey(ax.attributes,:subscenecontext)
        ax.attributes[:subscenecontext]=ctx
    end
    ax.attributes[:subscenecontext][]
end

"""
    setmakie!(MyMakie)

Set the Makie module.
"""
setmakie!(MyMakie) = global Makie=MyMakie



export multiscene, scene_interaction
export subscenecontext
export yieldwait
export add_scene!
export setmakie!
end
