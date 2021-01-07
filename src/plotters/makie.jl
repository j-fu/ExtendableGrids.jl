include("flippablelayout.jl")
using .FlippableLayout

function initialize_gridplot!(p::GridPlotContext,::Type{MakieType})
    Makie=p.context[:Plotter]
    if !isdefined(Makie.AbstractPlotting,:Box)
        error("Outdated version of AbstractPlotting. Please upgrade to at least 0.15.1")
    end
    FlippableLayout.setmakie!(Makie)
    layout=p.context[:layout]
    parent,flayout=flayoutscene(resolution=p.context[:resolution])
    p.context[:scene]=parent
    for I in CartesianIndices(layout)
        ctx=p.subplots[I]
        ctx[:figure]=parent
        ctx[:flayout]=flayout
    end
    p.context[:scene]
end

add_scene!(ctx,ax)=ctx[:flayout][ctx[:subplot]...]=ax


"""
     scene_interaction(update_scene,view,switchkeys::Vector{Symbol}=[:nothing])   

Control multiple scene elements via keyboard up/down keys. 
Each switchkey is assumed to correspond to one of these elements.
Pressing a switch key transfers control to its associated element.

Control of values of the current associated element is performed
by triggering change values via up/down (± 1)  resp. page_up/page_down (±10) keys

The update_scene callbac gets passed the change value and the symbol.
"""
function scene_interaction(update_scene,scene,Makie,switchkeys::Vector{Symbol}=[:nothing])
    function _inscene(scene,pos)
        area=scene.px_area[]
        pos[1]>area.origin[1] &&
            pos[1] < area.origin[1]+area.widths[1] &&
            pos[2]>area.origin[2] &&
            pos[2] < area.origin[2]+area.widths[2]
    end

    # Initial active switch key is the first in the vector passed
    activeswitch=Makie.Node(switchkeys[1])
    # Handle mouse position within scene-
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


function save(fname,p,::Type{MakieType})
    Makie=p.context[:Plotter]
    Makie.save(fname, p.context[:scene])
end


makestatus(grid::ExtendableGrid)="p: $(num_nodes(grid)) t: $(num_cells(grid)) b: $(num_bfaces(grid))"



#1D grid
function gridplot!(ctx, ::Type{MakieType}, ::Type{Val{1}}, grid)
    
    Makie=ctx[:Plotter]
    nregions=num_cellregions(grid)
    nbregions=num_bfaceregions(grid)
    
    function basemesh(grid)
        coord=vec(grid[Coordinates])
        xmin=minimum(coord)
        xmax=maximum(coord)
        h=(xmax-xmin)/40.0
        ncoord=length(coord)
        points=Vector{Point2f0}(undef,0)
        for i=1:ncoord
            push!(points,Point2f0(coord[i],h))
            push!(points,Point2f0(coord[i],-h))
        end
        points
    end

    function regionmesh(grid,iregion)
        coord=vec(grid[Coordinates])
        cn=grid[CellNodes]
        cr=grid[CellRegions]
        ncells=length(cr)
        points=Vector{Point2f0}(undef,0)
        for i=1:ncells
            if cr[i]==iregion
                push!(points,Point2f0(coord[cn[1,i]],0))
                push!(points,Point2f0(coord[cn[2,i]],0))
            end
        end
        points
    end

    function bmesh(grid,ibreg)
        coord=vec(grid[Coordinates])
        xmin=minimum(coord)
        xmax=maximum(coord)
        h=(xmax-xmin)/20.0
        nbfaces=num_bfaces(grid)
        bfacenodes=grid[BFaceNodes]
        bfaceregions=grid[BFaceRegions]
        points=Vector{Point2f0}(undef,0)
        for ibface=1:nbfaces
            if bfaceregions[ibface]==ibreg
                push!(points,Point2f0(coord[bfacenodes[1,ibface]],h))
                push!(points,Point2f0(coord[bfacenodes[1,ibface]],-h))
            end
        end
        points
    end

    
    if !haskey(ctx,:scene)
        ctx[:scene]=Makie.Axis(ctx[:figure],title=ctx[:title])
        ctx[:grid]=Makie.Node(grid)
        cmap=region_cmap(nregions)
        Makie.linesegments!(ctx[:scene],Makie.lift(g->basemesh(g), ctx[:grid]),color=:black)
        for i=1:nregions
            Makie.linesegments!(ctx[:scene],Makie.lift(g->regionmesh(g,i), ctx[:grid]) , color=cmap[i], strokecolor=:black,linewidth=4)
        end
        
        bcmap=bregion_cmap(nbregions)
        for i=1:nbregions
            Makie.linesegments!(ctx[:scene],Makie.lift(g->bmesh(g,i),ctx[:grid]), color=bcmap[i], linewidth=4)
        end
        add_scene!(ctx,ctx[:scene])
        Makie.display(ctx[:figure])
    else
        ctx[:grid][]=grid
        yieldwait()
    end
    ctx[:figure]
end

# 1D function
function gridplot!(ctx, ::Type{MakieType}, ::Type{Val{1}}, grid,func)
    Makie=ctx[:Plotter]

    if !haskey(ctx,:scene)
        ctx[:scene]=Makie.Axis(ctx[:figure],title=ctx[:title])
        ctx[:data]=Makie.Node((g=grid,f=func))
        Makie.lines!(ctx[:scene],Makie.lift(data->(vec(data.g[Coordinates]),data.f), ctx[:data]))
        add_scene!(ctx,ctx[:scene])
        Makie.display(ctx[:figure])
    else
        ctx[:data][]=(g=grid,f=func)
        yieldwait()
    end
    ctx[:figure]
end

# 2D grid
function gridplot!(ctx, ::Type{MakieType}, ::Type{Val{2}},grid)
    Makie=ctx[:Plotter]
    nregions=num_cellregions(grid)
    nbregions=num_bfaceregions(grid)


    
    if !haskey(ctx,:scene)
        ctx[:scene]=Makie.Axis(ctx[:figure],title=ctx[:title])
        ctx[:grid]=Makie.Node(grid)
        cmap=region_cmap(nregions)

        
        for i=1:nregions
            Makie.poly!(ctx[:scene],Makie.lift(g->regionmesh(g,i), ctx[:grid]) , color=cmap[i], strokecolor=:black)
        end

        bcmap=bregion_cmap(nbregions)
        for i=1:nbregions
            Makie.linesegments!(ctx[:scene],Makie.lift(g->bfacesegments(g,i),ctx[:grid]) , color=bcmap[i], linewidth=4)
        end
        add_scene!(ctx,ctx[:scene])
        Makie.display(ctx[:figure])
    else
        ctx[:grid][]=grid
        yieldwait()
    end
    ctx[:figure]
end



# 2D function
function gridplot!(ctx, ::Type{MakieType}, ::Type{Val{2}},grid, func)
    Makie=ctx[:Plotter]
    
    function make_mesh(grid::ExtendableGrid,func,elevation)
        coord=grid[Coordinates]
        npoints=num_nodes(grid)
        cellnodes=grid[CellNodes]
        if elevation ≈ 0.0
            points=[Point2f0(coord[1,i],coord[2,i]) for i=1:npoints]
        else
            points=[Point3f0(coord[1,i],coord[2,i],func[i]*elevation) for i=1:npoints]
        end
        faces=[TriangleFace(cellnodes[1,i],cellnodes[2,i],cellnodes[3,i]) for i=1:size(cellnodes,2)]
        Mesh(points,faces)
    end
    
    if !haskey(ctx,:scene)
        ctx[:scene]=Makie.Axis(ctx[:figure],title=ctx[:title])
        ctx[:data]=Makie.Node((g=grid,f=func,e=ctx[:elevation]))
        Makie.poly!(ctx[:scene],
                    Makie.lift(data->make_mesh(data.g,data.f,data.e), ctx[:data]),
                    color=Makie.lift(data->data.f,ctx[:data]),
                    colormap=ctx[:colormap])
        add_scene!(ctx,ctx[:scene])
        Makie.display(ctx[:figure])
    else
        ctx[:data][]=(g=grid,f=func,e=ctx[:elevation])
        yieldwait()
    end
    ctx[:figure]
end



# 3D Grid
function xyzminmax(grid::ExtendableGrid)
    coord=grid[Coordinates]
    ndim=size(coord,1)
    xyzmin=zeros(ndim)
    xyzmax=ones(ndim)
    for idim=1:ndim
        @views mn,mx=extrema(coord[idim,:])
        xyzmin[idim]=mn
        xyzmax[idim]=mx
    end
    xyzmin,xyzmax
end

const keyboardhelp=
"""
Keyboard interactions:
          x: control xplane
          y: control yplane
          z: control zplane
          l: control isolevel
    up/down: fine control control value
pgup/pgdown: coarse control control value
          h: print this message
"""
function gridplot!(ctx, ::Type{MakieType}, ::Type{Val{3}}, grid)

    make_mesh(pts,fcs)=Mesh(meta(pts,normals=normals(pts, fcs)),fcs)
    
    nregions=num_cellregions(grid)
    nbregions=num_bfaceregions(grid)
    
    Makie=ctx[:Plotter]
    xyzmin,xyzmax=xyzminmax(grid)
    xyzstep=(xyzmax-xyzmin)/100
    
    function adjust_planes()
        ctx[:xplane]=max(xyzmin[1],min(xyzmax[1],ctx[:xplane]) )
        ctx[:yplane]=max(xyzmin[2],min(xyzmax[2],ctx[:yplane]) )
        ctx[:zplane]=max(xyzmin[3],min(xyzmax[3],ctx[:zplane]) )
    end

    adjust_planes()

    

    if !haskey(ctx,:scene)
        ctx[:scene]=Makie.LScene(ctx[:figure])
        ctx[:data]=Makie.Node((g=grid,x=ctx[:xplane],y=ctx[:yplane],z=ctx[:zplane]))
        ctx[:facedata]=Makie.lift(
            d->extract_visible_bfaces3D(d.g,
                                        (d.x,d.y,d.z),
                                        primepoints=hcat(xyzmin,xyzmax),
                                        Tp=Point3f0,
                                        Tf=GLTriangleFace),
            ctx[:data])
        
        ctx[:facemeshes]=Makie.lift(d->[make_mesh(d[1][i],d[2][i]) for i=1:nbregions], ctx[:facedata])
        
        if ctx[:interior]
            cmap=region_cmap(nregions)
            ctx[:celldata]=Makie.lift(
                d->extract_visible_cells3D(d.g,
                                           (d.x,d.y,d.z),
                                           primepoints=hcat(xyzmin,xyzmax),
                                           Tp=Point3f0,
                                           Tf=GLTriangleFace),
                ctx[:data])
            ctx[:cellmeshes]=Makie.lift(d->[make_mesh(d[1][i],d[2][i]) for i=1:nregions], ctx[:celldata])
            for i=1:nregions
                Makie.mesh!(ctx[:scene],Makie.lift(d->d[i], ctx[:cellmeshes]),
                            color=cmap[i],
                            backlight=1f0
                            )
                if ctx[:edges]
                    Makie.wireframe!(ctx[:scene],Makie.lift(d->d[i], ctx[:cellmeshes]),strokecolor=:black)
                end
            end
        end

        bcmap=bregion_cmap(nbregions)
        for i=1:nbregions
            Makie.mesh!(ctx[:scene],Makie.lift(d->d[i], ctx[:facemeshes]),
                        color=bcmap[i],
                        backlight=1f0
                        )
            if ctx[:edges]
                Makie.wireframe!(ctx[:scene],Makie.lift(d->d[i], ctx[:facemeshes]),strokecolor=:black)
            end
        end
        
        
        scene_interaction(ctx[:scene].scene,Makie,[:z,:y,:x]) do delta,key
            if key==:x
                ctx[:xplane]+=delta*xyzstep[1]
            elseif key==:y
                ctx[:yplane]+=delta*xyzstep[2]
            elseif key==:z
                ctx[:zplane]+=delta*xyzstep[3]
            end
            adjust_planes()
            ctx[:data][]=(g=grid,x=ctx[:xplane],y=ctx[:yplane],z=ctx[:zplane])
        end
        add_scene!(ctx,ctx[:scene])
        Makie.display(ctx[:figure])
    else
        ctx[:data][]=(g=grid,x=ctx[:xplane],y=ctx[:yplane],z=ctx[:zplane])
        yieldwait()
    end
    ctx[:figure]
end

function gridplot!(ctx, ::Type{MakieType}, ::Type{Val{3}}, grid , func)
    
    make_mesh(pts,fcs)=Mesh(pts,fcs)
    
    function make_mesh(pts,fcs,vals)
        colors = Makie.AbstractPlotting.interpolated_getindex.((cmap,), vals, (fminmax,))
        GeometryBasics.Mesh(meta(pts, color=colors,normals=normals(pts, fcs)), fcs)
    end
        
    
    nregions=num_cellregions(grid)
    nbregions=num_bfaceregions(grid)
    
    Makie=ctx[:Plotter]
    cmap = Makie.to_colormap(ctx[:colormap])
    xyzmin,xyzmax=xyzminmax(grid)
    xyzstep=(xyzmax-xyzmin)/100

    fminmax=extrema(func)
    fstep=(fminmax[2]-fminmax[1])/100
    if fstep≈0
        fstep=0.1
    end
    
    function adjust_planes()
        ctx[:xplane]=max(xyzmin[1],min(xyzmax[1],ctx[:xplane]) )
        ctx[:yplane]=max(xyzmin[2],min(xyzmax[2],ctx[:yplane]) )
        ctx[:zplane]=max(xyzmin[3],min(xyzmax[3],ctx[:zplane]) )
        ctx[:flevel]=max(fminmax[1],min(fminmax[2],ctx[:flevel]))
    end
    
    adjust_planes()


    makeplanes(x,y,z)=[[1,0,0,-x], 
                       [0,1,0,-y], 
                       [0,0,1,-z]]
    
    if !haskey(ctx,:scene)
        ctx[:scene]=Makie.LScene(ctx[:figure])
        ctx[:data]=Makie.Node((g=grid,f=func,x=ctx[:xplane],y=ctx[:yplane],z=ctx[:zplane],l=ctx[:flevel]))
                                  
        if ctx[:outline]


            ctx[:facedata]=Makie.lift(d->extract_visible_bfaces3D(d.g,
                                                                  xyzmax,
                                                                  primepoints=hcat(xyzmin,xyzmax),
                                                                  Tp=Point3f0,
                                                                  Tf=GLTriangleFace),
                                      ctx[:data])
            ctx[:facemeshes]=Makie.lift(d->[make_mesh(d[1][i],d[2][i]) for i=1:nbregions], ctx[:facedata])
            bcmap=bregion_cmap(nbregions)
            for i=1:nbregions
                Makie.mesh!(ctx[:scene],Makie.lift(d->d[i], ctx[:facemeshes]),
                color=(bcmap[i],ctx[:alpha]),
                transparency=true,
                backlight=1f0
                )
            end
        end
        
        Makie.mesh!(ctx[:scene],
                    Makie.lift(d->make_mesh(marching_tetrahedra(d.g,
                                                                d.f,
                                                                makeplanes(d.x,d.y,d.z),
                                                                [d.l],
                                                                primepoints=hcat(xyzmin,xyzmax),
                                                                primevalues=fminmax,
                                                                Tp=Point3f0,
                                                                Tf=GLTriangleFace,
                                                                Tv=Float32)...),ctx[:data]),
                    backlight=1f0)
        add_scene!(ctx,ctx[:scene])
        Makie.display(ctx[:figure])
        
        scene_interaction(ctx[:scene].scene,Makie,[:z,:y,:x,:l]) do delta,key
            if key==:x
                ctx[:xplane]+=delta*xyzstep[1]
            elseif key==:y
                ctx[:yplane]+=delta*xyzstep[2]
            elseif key==:z
                ctx[:zplane]+=delta*xyzstep[3]
            elseif key==:l
                ctx[:flevel]+=delta*fstep
            end
            adjust_planes()
            ctx[:data][]=(g=grid,f=func,x=ctx[:xplane],y=ctx[:yplane],z=ctx[:zplane],l=ctx[:flevel])
        end
    else
        ctx[:data][]=(g=grid,f=func,x=ctx[:xplane],y=ctx[:yplane],z=ctx[:zplane],l=ctx[:flevel])
        yieldwait()
    end
    ctx[:figure]
end






        # TODO: allow aspect scaling
        # if ctx[:aspect]>1.0
        #     Makie.scale!(ctx[:scene],ctx[:aspect],1.0)
        # else
        #     Makie.scale!(ctx[:scene],1.0,1.0/ctx[:aspect])
        # end

        # TODO: use distinguishable colors
        # http://juliagraphics.github.io/Colors.jl/stable/colormapsandcolorscales/#Generating-distinguishable-colors-1
    
# TODO: a priori angles aka pyplot3D
# rect = ctx[:scene]
# azim=ctx[:azim]
# elev=ctx[:elev]
# arr = normalize([cosd(azim/2), 0, sind(azim/2), -sind(azim/2)])
# Makie.rotate!(rect, Makie.Quaternionf0(arr...))

