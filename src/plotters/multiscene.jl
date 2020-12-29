module MultiScene
# (c) Julius Krummbiegel, Jürgen Fuhrmann


Makie=nothing

setmakie(MyMakie) = global Makie=MyMakie



###############################################################
#
# Generic part

"""
    Check if position is  within pixel area of scene
"""
function inscene(scene,pos)
    area=scene.px_area[]
    pos[1]>area.origin[1] &&
        pos[1] < area.origin[1]+area.widths[1] &&
        pos[2]>area.origin[2] &&
        pos[2] < area.origin[2]+area.widths[2]
end



"""
Control multiple scene elements via keyboard up/down keys. 
Each switchkey is assumed to correspond to one of these elements.
Pressing a switch key transfers control via up/down resp. page_up/page_down
to its associated element.
"""
function scene_interaction(update_scene,scene,switchkeys::Vector{Symbol}=[:nothing])
    mouseposition=Makie.Node((0,0))
    activeswitch=Makie.Node(switchkeys[1])
    
    Makie.on(m->mouseposition[]=m, scene.events.mouseposition)

    Makie.on(scene.events.keyboardbuttons) do buttons
        if inscene(scene,mouseposition[])
            for i=1:length(switchkeys)
                if switchkeys[i]!=:nothing&&Makie.ispressed(scene,getproperty(Makie.Keyboard,switchkeys[i]))
                    activeswitch[]=switchkeys[i]
                    update_scene(0,switchkeys[i])
                    return
                end
            end
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
Combine scene with title header
"""
textposition(scene)=Makie.lift(a->Makie.Vec2f0(Makie.widths(a)[1] ./ 2 , 0), Makie.pixelarea(scene))
header(scene,txt)=Makie.text(txt,textsize = 20,raw=true,position=textposition(scene),camera=Makie.campixel!,align = (:center, :bottom))
footer(scene,status)=Makie.text(Makie.lift(a->a,status),textsize = 15,raw=true,position=textposition(scene),camera=Makie.campixel!,align = (:center, :bottom))
function add_scene!(ax,scene;title=" ",status=nothing)
    if isnothing(status)
        Makie.hbox(scene,header(scene,title),parent=ax.scene)
    else
        Makie.hbox(footer(scene,status),scene,header(scene,title),parent=ax.scene)
    end
end


#####################################################
# Handle blocking/unblocking via space key
mutable struct Blocker
    condition::Condition
    blocked::Bool
    Blocker(;blocked=false)=new(Condition(),blocked)
end


# block/unblock, and notify about unblocking
function toggle(blocker::Blocker)
    if blocker.blocked
        blocker.blocked=false
        notify(blocker.condition)
    else
        blocker.blocked=true
    end
end

#
# If blocked, wait for notification.
# The yield goes here conveniently as well.
#
function waitblocker(blocker::Blocker)
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

waitblocker()=waitblocker(globalblocker)


"""
Create a scene with given layout grid. Returns an array of subscenes
according to the layout parameter.

The `,` key switches between focused view showing only one subscene
and "gallery view" showing all subscenes at once.
"""
function multiscene(;layout=(1,1), blocked=false, resolution=(500,500))
    (scene, makielayout) = Makie.layoutscene(resolution = resolution)
    offscreen_gl = Makie.GridLayout(bbox = Makie.BBox(-500, -400, -500, -400))
    axs = makielayout[] = [Makie.LScene(scene) for _ in CartesianIndices(layout)]

    gallery_view=Makie.Node(true)
    mouseposition=Makie.Node((0,0))
    globalblocker.blocked=blocked
    
    Makie.on(m->mouseposition[]=m, scene.events.mouseposition)
    
    function focus(i)
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
    
    function show_all()
        for idx in CartesianIndices(axs)
            makielayout[Tuple(idx)...]=  axs[Tuple(idx)...]
        end
        gallery_view[]=true
    end
    
    function child(mouseposition)
        for i=1:length(scene.children)
            if inscene(scene.children[i],mouseposition)
                return i
            end
        end
        return 0
    end
    
    
    Makie.on(scene.events.keyboardbuttons) do buttons
        if Makie.ispressed(scene, Makie.Keyboard.comma)
            gallery_view[] ? focus(child(mouseposition[])) : show_all() 
        end
        if Makie.ispressed(scene, Makie.Keyboard.space)
            toggle(globalblocker)
        end
    end
    scene,axs
end

function subscenecontext(ax; ctx=Dict{Symbol,Any}())
    if !haskey(ax.attributes,:subscenecontext)
        ax.attributes[:subscenecontext]=ctx
    end
    ax.attributes[:subscenecontext][]
end

export multiscene, scene_interaction
export subscenecontext,waitblocker
export add_scene!
end # Module MultiScene2
