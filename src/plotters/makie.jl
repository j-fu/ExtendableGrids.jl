
function initialize_context!(ctx::PlotterContext,::Type{MakieType})
    ctx
end

function make_scene!(ctx::PlotterContext)
    Makie=ctx[:Plotter]
    if !haskey(ctx,:scene)
       ctx[:scene]=Makie.Scene(scale_plot=false, resolution=ctx[:resolution])
    end
    ctx
end


points(coord) =points(coord,Val{size(coord,1)})
points(coord,::Type{Val{2}}) =[Point2f0(coord[1,i],coord[2,i]) for i=1:size(coord,2)]
points(coord,::Type{Val{3}}) =[Point3f0(coord[1,i],coord[2,i],coord[3,i]) for i=1:size(coord,2)]

faces(cellnodes) = faces(cellnodes,Val{size(cellnodes,1)-1})
faces(cellnodes,::Type{Val{2}}) = [TriangleFace(cellnodes[1,i],cellnodes[2,i],cellnodes[3,i]) for i=1:size(cellnodes,2)]

make_mesh(grid::ExtendableGrid)=Mesh(points(grid[Coordinates]),faces(grid[CellNodes]))



function make_mesh(grid::ExtendableGrid, elevation; elevation_factor=1.0)
    coord=grid[Coordinates]
    npoints=num_nodes(grid)
    points=[Point3f0(coord[1,i],coord[2,i],elevation[i]*elevation_factor) for i=1:npoints]
    Mesh(points,faces(grid[CellNodes]))
end

         
function region_bfacesegments(grid::ExtendableGrid,ibreg)
    coord=grid[Coordinates]
    nbfaces=num_bfaces(grid)
    bfacenodes=grid[BFaceNodes]
    bfaceregions=grid[BFaceRegions]
    points=Vector{Point2f0}(undef,0)
    for ibface=1:nbfaces
        if bfaceregions[ibface]==ibreg
            push!(points,Point2f0(coord[1,bfacenodes[1,ibface]],coord[2,bfacenodes[1,ibface]]))
            push!(points,Point2f0(coord[1,bfacenodes[2,ibface]],coord[2,bfacenodes[2,ibface]]))
        end
    end
    points
end


function plot!(ctx, ::Type{MakieType}, ::Type{Val{2}},grid)
    Makie=ctx[:Plotter]
    nregions=num_cellregions(grid)
    if nregions==1
        meshes=[make_mesh(grid)]
    else
        meshes=[make_mesh(subgrid(grid,[iregion])) for iregion=1:nregions]
    end

    make_scene!(ctx)
    nbregions=num_bfaceregions(grid)
    bsegments=[region_bfacesegments(grid,iregion) for iregion=1:nbregions]
    
    if !haskey(ctx,:meshes)|| ctx[:num_cellregions]!=nregions
        ctx[:num_cellregions]=nregions
        # if ctx[:aspect]>1.0
        #     Makie.scale!(ctx[:scene],ctx[:aspect],1.0)
        # else
        #     Makie.scale!(ctx[:scene],1.0,1.0/ctx[:aspect])
        # end
        ctx[:meshes]=[Makie.Node(meshes[i]) for i=1:nregions]
        ctx[:bsegments]=[Makie.Node(bsegments[i]) for i=1:nbregions]
        for i=1:nregions
            Makie.poly!(ctx[:scene],Makie.lift(a->a, ctx[:meshes][i]) , color=frgb(Makie,i,nregions,pastel=true), strokecolor=:black)
        end
        for i=1:nbregions
            Makie.linesegments!(ctx[:scene],Makie.lift(a->a, ctx[:bsegments][i]) , color=frgb(Makie,i,nbregions), linewidth=4)
        end
        Makie.display(ctx[:scene])
    else
        for i=1:nregions
            ctx[:meshes][i][]=meshes[i]
        end
        for i=1:nbregions
            ctx[:bsegments][i][]=bsegments[i]
        end
        if ctx[:show]
            Makie.update!(ctx[:scene])
        end
        Makie.sleep(1.0e-10)
    end
    ctx[:scene]
end


function plot!(ctx, ::Type{MakieType}, ::Type{Val{2}},grid, func)
    Makie=ctx[:Plotter]
    make_scene!(ctx)
    if ctx[:elevation]!=0
        mesh=make_mesh(grid,func,elevation_factor=ctx[:elevation])
    else
        mesh=make_mesh(grid)
    end        
    if !haskey(ctx,:meshnode)
        ctx[:meshnode]=Makie.Node(mesh)
        ctx[:colornode]=Makie.Node(func)
        Makie.poly!(ctx[:scene],Makie.lift(a->a, ctx[:meshnode]) , color=Makie.lift(a->a, ctx[:colornode]), colormap=ctx[:colormap])
        Makie.display(ctx[:scene])
    else
        ctx[:meshnode][]=mesh
        ctx[:colornode][]=func
        if ctx[:show]
            Makie.update!(ctx[:scene])
        end
        Makie.sleep(1.0e-2)
    end
    ctx[:scene]
end


#####################################################################################
# 3D
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


function set_interaction_handler(ctx,update_scene,xyzmin,xyzstep,xyzmax)
    Makie=ctx[:Plotter]
    ctx[:moving]=:zplane
    ctx[:movdir]=3
    # used later when we have layouting
    on(ctx[:scene].events.mouseposition) do mpos
        ctx[:mouseposition][]=mpos
    end
    
    on(ctx[:scene].events.keyboardbuttons) do buttons
        scene= ctx[:scene]
        mpos=ctx[:mouseposition][]
        # area=ctx[:scene].px_area
        # if (mpos[1]>area.origin[1] && mpos[1] < area.origin[1]+area.widths[1]
        #     && mpos[2]>area.origin[2] && mpos[2] < area.origin[2]+area.widths[2]
        #     )
        if Makie.ispressed(scene, Makie.Keyboard.x)
            ctx[:moving]=:xplane
            ctx[:movdir]=1
        elseif Makie.ispressed(scene, Makie.Keyboard.y)
            ctx[:moving]=:yplane
            ctx[:movdir]=2
        elseif Makie.ispressed(scene, Makie.Keyboard.z)
            ctx[:moving]=:zplane
            ctx[:movdir]=3
        elseif Makie.ispressed(scene, Makie.Keyboard.l)
            ctx[:moving]=:flevel
            ctx[:movdir]=4
        elseif Makie.ispressed(scene, Makie.Keyboard.down)
            ctx[ctx[:moving]]-=xyzstep[ctx[:movdir]]
            ctx[ctx[:moving]]=max(xyzmin[ctx[:movdir]],ctx[ctx[:moving]])
            ctx[:footer][]=@sprintf("%s=%.2e",ctx[:moving],ctx[ctx[:moving]])
            update_scene([ctx[:xplane],ctx[:yplane],ctx[:zplane],ctx[:flevel]])
        elseif Makie.ispressed(scene, Makie.Keyboard.page_down)
            ctx[ctx[:moving]]-=xyzstep[ctx[:movdir]]*10
            ctx[ctx[:moving]]=max(xyzmin[ctx[:movdir]],ctx[ctx[:moving]])
            ctx[:footer][]=@sprintf("%s=%.2e",ctx[:moving],ctx[ctx[:moving]])
            update_scene([ctx[:xplane],ctx[:yplane],ctx[:zplane],ctx[:flevel]])
        elseif Makie.ispressed(scene, Makie.Keyboard.up)
            ctx[ctx[:moving]]+=xyzstep[ctx[:movdir]]
            ctx[ctx[:moving]]=min(xyzmax[ctx[:movdir]],ctx[ctx[:moving]])
            ctx[:footer][]=@sprintf("%s=%.2e",ctx[:moving],ctx[ctx[:moving]])
            update_scene([ctx[:xplane],ctx[:yplane],ctx[:zplane],ctx[:flevel]])
        elseif Makie.ispressed(scene, Makie.Keyboard.page_up)
            ctx[ctx[:moving]]+=xyzstep[ctx[:movdir]]*10
            ctx[ctx[:moving]]=min(xyzmax[ctx[:movdir]],ctx[ctx[:moving]])
            ctx[:footer][]=@sprintf("%s=%.2e",ctx[:moving],ctx[ctx[:moving]])
            update_scene([ctx[:xplane],ctx[:yplane],ctx[:zplane],ctx[:flevel]])
        elseif Makie.ispressed(scene, Makie.Keyboard.h)
            println(keyboardhelp)
        end
    end
end



function plot!(ctx, ::Type{MakieType}, ::Type{Val{3}},grid)
    
    function make_mesh(pts,fcs)
        pts=points(pts)
        # reinterpret_array kann kein push... so we need to seed it befor
        # so that it is not empty
        push!(pts,Point3f0(xyzmin...))
        push!(pts,Point3f0(xyzmax...))
        fcs=faces(fcs)
        Mesh(meta(pts,normals=normals(pts, fcs)),fcs)
    end

    
    
    Makie=ctx[:Plotter]
    ctx[:mouseposition]=Makie.Node((0,0))

    xyzmin,xyzmax=xyzminmax(grid)

    xyzstep=(xyzmax-xyzmin)/100
    
    ctx[:xplane]=min(xyzmax[1],ctx[:xplane])
    ctx[:yplane]=min(xyzmax[2],ctx[:yplane])
    ctx[:zplane]=min(xyzmax[3],ctx[:zplane])

    ctx[:xplane]=max(xyzmin[1],ctx[:xplane])
    ctx[:yplane]=max(xyzmin[2],ctx[:yplane])
    ctx[:zplane]=max(xyzmin[3],ctx[:zplane])

    
    function create_scene(xyzcut)
        alpha=ctx[:alpha]
        trans=alpha<1
        nregions=num_cellregions(grid)
        nbregions=num_bfaceregions(grid)
        coord=grid[Coordinates]
        # TODO: allow aspect scaling
        # if ctx[:aspect]>1.0
        #     Makie.scale!(ctx[:scene],ctx[:aspect],1.0)
        # else
        #     Makie.scale!(ctx[:scene],1.0,1.0/ctx[:aspect])
        # end
        
        scene = Makie.Scene(scale_plot=false)
    
        regpoints,regfacets=extract_visible_cells3D(grid,xyzcut)
        bregpoints,bregfacets=extract_visible_bfaces3D(grid,xyzcut)
        ctx[:nregions]=nregions
        ctx[:nbregions]=nbregions
        ctx[:meshes]=   [Makie.Node(make_mesh(regpoints[iregion],regfacets[iregion])) for iregion=1:nregions]
        ctx[:bsegments]=[Makie.Node(make_mesh(bregpoints[ibregion],bregfacets[ibregion])) for ibregion=1:nbregions]
        
        # TODO: use distinguishable colors
        # http://juliagraphics.github.io/Colors.jl/stable/colormapsandcolorscales/#Generating-distinguishable-colors-1
        if ctx[:interior]
            for i=1:nregions
                Makie.mesh!(scene,Makie.lift(a->a, ctx[:meshes][i]),
                            color=(frgb(Makie,i,nregions+nbregions,pastel=true),alpha),
                            backlight=1f0,
                            transparency=trans
                            )
                if (!trans)
                    Makie.wireframe!(scene,Makie.lift(a->a, ctx[:meshes][i]),strokecolor=:black)
                end
            end
        end
        for i=1:nbregions
            Makie.mesh!(scene,Makie.lift(a->a, ctx[:bsegments][i]),
                        color=(frgb(Makie,nregions+i,nregions+nbregions),alpha),
                        transparency=trans)
            if (!trans)
                Makie.wireframe!(scene,Makie.lift(a->a, ctx[:bsegments][i]) , strokecolor=:black)
            end
        end


        pos = Makie.lift(Makie.pixelarea(scene)) do area
            x = widths(area)[1] ./ 2
            Vec2f0(x,0) # offset 10px, to give it some space
        end

        ctx[:footer]=Makie.Node(" ")
        
        title = Makie.text(ctx[:title],textsize = 20,raw=true,position=pos,camera=Makie.campixel!,align = (:center, :bottom))
        footer = Makie.text(Makie.lift(a->a,ctx[:footer]),textsize = 15,raw=true,position=pos,camera=Makie.campixel!,align = (:center, :bottom))


        ctx[:scene]=Makie.hbox(footer,scene,title)
              

        set_interaction_handler(ctx,update_scene,xyzmin,xyzstep,xyzmax)
        
        Makie.display(ctx[:scene])
    end

    function update_scene(xyzcut)
        if ctx[:interior]
            regpoints,regfacets=extract_visible_cells3D(grid,xyzcut)
            for i=1:ctx[:nregions]
                ctx[:meshes][i][]=make_mesh(regpoints[i],regfacets[i])
            end
        end
        bregpoints,bregfacets=extract_visible_bfaces3D(grid,xyzcut)
        for i=1:ctx[:nbregions]
            ctx[:bsegments][i][]=make_mesh(bregpoints[i],bregfacets[i])
        end
    end
    
    xyz=[ctx[:xplane],ctx[:yplane],ctx[:zplane]]
    
    if (!haskey(ctx,:scene)||
        ctx[:nregions]!=nregions ||
        ctx[:nbregions]!=nbregions)
        create_scene(xyz)
    else
        update_scene(xyz)
    end
    ctx[:scene]
end



function plot!(ctx, ::Type{MakieType}, ::Type{Val{3}}, grid , func)

    function make_mesh(pts,fcs)
        pts=points(pts)
        push!(pts,Point3f0(xyzmin[1:3]...))
        push!(pts,Point3f0(xyzmax[1:3]...))
        fcs=faces(fcs)
        Mesh(pts,fcs)
    end
    
    
    
    function make_mesh(pts,fcs,vals)
        pts=points(pts)
        push!(pts,Point3f0(xyzmin[1:3]...))
        push!(pts,Point3f0(xyzmax[1:3]...))
        fcs=faces(fcs)
        push!(vals,0)
        push!(vals,0)
        colors = Makie.AbstractPlotting.interpolated_getindex.((cmap,), vals, (fminmax,))
        GeometryBasics.Mesh(meta(pts, color=colors,normals=normals(pts, fcs)), fcs)
    end


    Makie=ctx[:Plotter]
    ctx[:mouseposition]=Makie.Node((0,0))

    xyzmin,xyzmax=xyzminmax(grid)

    xyzstep=(xyzmax-xyzmin)/100

    fminmax=extrema(func)
    fstep=(fminmax[2]-fminmax[1])/100
    if fstep≈0
        fstep=0.1
    end
    xyzmin=[xyzmin...,fminmax[1]]
    xyzmax=[xyzmax...,fminmax[2]]
    xyzstep=[xyzstep...,fstep]
    
    
    cmap = Makie.to_colormap(ctx[:colormap])

    ctx[:xplane]=min(xyzmax[1],ctx[:xplane])
    ctx[:yplane]=min(xyzmax[2],ctx[:yplane])
    ctx[:zplane]=min(xyzmax[3],ctx[:zplane])
    ctx[:flevel]=min(xyzmax[4],ctx[:flevel])

    ctx[:xplane]=max(xyzmin[1],ctx[:xplane])
    ctx[:yplane]=max(xyzmin[2],ctx[:yplane])
    ctx[:zplane]=max(xyzmin[3],ctx[:zplane])
    ctx[:flevel]=max(xyzmin[4],ctx[:flevel])

    
    ctx[:grid]=grid
    ctx[:func]=func

    makeplanes(xyzf)=[ [1,0,0,-xyzf[1]], 
                       [0,1,0,-xyzf[2]], 
                       [0,0,1,-xyzf[3]]]
    
    function create_scene(xyzf)
        alpha=ctx[:alpha]
        nbregions=num_bfaceregions(grid)
        nregions=num_cellregions(grid)

        coord=grid[Coordinates]
            
        scene = Makie.Scene(scale_plot=false)
        ctx[:num_cellregions]=nregions
        ctx[:num_bfaceregions]=nbregions
        
        bregpoints,bregfacets=extract_visible_bfaces3D(grid,xyzmax)
        if ctx[:outline]
            ctx[:bsegments]=[Makie.Node(make_mesh(bregpoints[ibregion],bregfacets[ibregion])) for ibregion=1:nbregions]
        end
        
        ppoints,pfacets,pvalues=marching_tetrahedra(grid,func,makeplanes(xyzf),[xyzf[4]])
        ctx[:pmesh]=Makie.Node(make_mesh(ppoints,pfacets,pvalues))
        Makie.mesh!(scene, Makie.lift(a->a, ctx[:pmesh]),backlight=1f0)

        if ctx[:outline]
            for i=1:ctx[:num_bfaceregions]
                Makie.mesh!(scene,Makie.lift(a->a, ctx[:bsegments][i]),
                            color=(frgb(Makie,nregions+i,nregions+nbregions),ctx[:alpha]),
                            transparency=true)
            end
        end

        pos = Makie.lift(Makie.pixelarea(scene)) do area
            x = widths(area)[1] ./ 2
            Vec2f0(x,0) # offset 10px, to give it some space
        end

        ctx[:footer]=Makie.Node(" ")
        
        title = Makie.text(ctx[:title],textsize = 20,raw=true,position=pos,camera=Makie.campixel!,align = (:center, :bottom))
        footer = Makie.text(Makie.lift(a->a,ctx[:footer]),textsize = 15,raw=true,position=pos,camera=Makie.campixel!,align = (:center, :bottom))
        
        
        ctx[:scene]=Makie.hbox(footer,scene,title)
        
        
        set_interaction_handler(ctx,update_scene,xyzmin,xyzstep,xyzmax)
        Makie.display(ctx[:scene])

    end

    
    function update_scene(xyzf)
        ppoints,pfacets,pvalues=marching_tetrahedra(grid,func,makeplanes(xyzf),[xyzf[4]])
        ctx[:pmesh][]= make_mesh(ppoints,pfacets,pvalues)
    end

    xyzf=[ctx[:xplane],ctx[:yplane],ctx[:zplane],ctx[:flevel]]
    
    if (!haskey(ctx,:scene)||ctx[:num_bfaceregions]!=nbregions )
        create_scene(xyzf)
    else
        update_scene(xyzf)
    end
    ctx[:scene]
end



# TODO: allow aspect scaling
# if ctx[:aspect]>1.0
#     Makie.scale!(ctx[:scene],ctx[:aspect],1.0)
# else
#     Makie.scale!(ctx[:scene],1.0,1.0/ctx[:aspect])
# end

# TODO: a priori angles aka pyplot3D
# rect = ctx[:scene]
# azim=ctx[:azim]
# elev=ctx[:elev]
# arr = normalize([cosd(azim/2), 0, sind(azim/2), -sind(azim/2)])
# Makie.rotate!(rect, Makie.Quaternionf0(arr...))
