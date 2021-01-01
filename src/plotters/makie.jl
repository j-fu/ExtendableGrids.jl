include("multiscene.jl")
using .MultiScene

function initialize_plot!(p::PlotContext,::Type{MakieType})
    Makie=p.context[:Plotter]
    MultiScene.setmakie!(Makie)
    layout=p.context[:layout]
    parent,subscenes=multiscene(layout=layout,resolution=p.context[:resolution])
    p.context[:scene]=parent
    for I in CartesianIndices(layout)
        ctx=p.subplots[I]
        ctx[:ax]=subscenes[I]
        ctx[:figure]=parent
    end
    p.context[:scene]
end


makestatus(grid::ExtendableGrid)="p: $(num_nodes(grid)) t: $(num_cells(grid)) b: $(num_bfaces(grid))"



#1D grid
function plot!(ctx, ::Type{MakieType}, ::Type{Val{1}}, grid)
    
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
        ctx[:scene]=Makie.Scene(scale_plot=false)
        ctx[:grid]=Makie.Node(grid)
        Makie.linesegments!(ctx[:scene],Makie.lift(g->basemesh(g), ctx[:grid]),color=:black)
        for i=1:nregions
            Makie.linesegments!(ctx[:scene],Makie.lift(g->regionmesh(g,i), ctx[:grid]) , color=frgb(Makie,i,nregions,pastel=true), strokecolor=:black)
        end
        
        for i=1:nbregions
            Makie.linesegments!(ctx[:scene],Makie.lift(g->bmesh(g,i),ctx[:grid]) , color=frgb(Makie,i,nbregions), linewidth=4)
        end
        add_scene!(ctx[:ax],ctx[:scene],title=ctx[:title],status=Makie.lift(g->makestatus(g),ctx[:grid]))
        Makie.display(ctx[:figure])
    else
        ctx[:grid][]=grid
        yieldwait()
    end
    ctx[:figure]
end

# 1D function
function plot!(ctx, ::Type{MakieType}, ::Type{Val{1}}, grid,func)
    Makie=ctx[:Plotter]

    if !haskey(ctx,:scene)
        ctx[:scene]=Makie.Scene(scale_plot=true)
        ctx[:data]=Makie.Node((g=grid,f=func))
        Makie.lines!(ctx[:scene],Makie.lift(data->(vec(data.g[Coordinates]),data.f), ctx[:data]))

        add_scene!(ctx[:ax],ctx[:scene],title=ctx[:title])
        Makie.display(ctx[:figure])
    else
        ctx[:data][]=(g=grid,f=func)
        yieldwait()
    end
    ctx[:figure]
end

# 2D grid
function plot!(ctx, ::Type{MakieType}, ::Type{Val{2}},grid)
    Makie=ctx[:Plotter]
    nregions=num_cellregions(grid)
    nbregions=num_bfaceregions(grid)

    function regionmesh(grid,iregion)
        coord=grid[Coordinates]
        cn=grid[CellNodes]
        cr=grid[CellRegions]
        @views points=[Point2f0(coord[:,i]) for i=1:size(coord,2)]
        faces=Vector{GLTriangleFace}(undef,0)
        for i=1:length(cr)
            if cr[i]==iregion
                @views push!(faces,cn[:,i])
            end
        end
        Mesh(points,faces)
    end

    function bfacesegments(grid,ibreg)
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

    
    if !haskey(ctx,:scene)
        ctx[:scene]=Makie.Scene(scale_plot=false)
        ctx[:grid]=Makie.Node(grid)
        for i=1:nregions
            Makie.poly!(ctx[:scene],Makie.lift(g->regionmesh(g,i), ctx[:grid]) , color=frgb(Makie,i,nregions,pastel=true), strokecolor=:black)
        end
        for i=1:nbregions
            Makie.linesegments!(ctx[:scene],Makie.lift(g->bfacesegments(g,i),ctx[:grid]) , color=frgb(Makie,i,nbregions), linewidth=4)
        end
        add_scene!(ctx[:ax],ctx[:scene],title=ctx[:title],status=Makie.lift(g->makestatus(g),ctx[:grid]))
        Makie.display(ctx[:figure])
    else
        ctx[:grid][]=grid
        yieldwait()
    end
    ctx[:figure]
end



# 2D function
function plot!(ctx, ::Type{MakieType}, ::Type{Val{2}},grid, func)
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
        ctx[:scene]=Makie.Scene(scale_plot=false)
        ctx[:data]=Makie.Node((g=grid,f=func,e=ctx[:elevation]))
        Makie.poly!(ctx[:scene],
                    Makie.lift(data->make_mesh(data.g,data.f,data.e), ctx[:data]),
                    color=Makie.lift(data->data.f,ctx[:data]),
                    colormap=ctx[:colormap])

        add_scene!(ctx[:ax],ctx[:scene],title=ctx[:title])
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
function plot!(ctx, ::Type{MakieType}, ::Type{Val{3}}, grid)

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
        ctx[:scene]=Makie.Scene(scale_plot=false)
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
                            color=(frgb(Makie,i,nregions,pastel=true)),
                            backlight=1f0
                            )
                if ctx[:edges]
                    Makie.wireframe!(ctx[:scene],Makie.lift(d->d[i], ctx[:cellmeshes]),strokecolor=:black)
                end
            end
        end

        for i=1:nbregions
            Makie.mesh!(ctx[:scene],Makie.lift(d->d[i], ctx[:facemeshes]),
                        color=(frgb(Makie,i,nbregions,pastel=false)),
                        backlight=1f0
                        )
            if ctx[:edges]
                Makie.wireframe!(ctx[:scene],Makie.lift(d->d[i], ctx[:facemeshes]),strokecolor=:black)
            end
        end
        
        
        scene_interaction(ctx[:scene],[:z,:y,:x]) do delta,key
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
        add_scene!(ctx[:ax],ctx[:scene],title=ctx[:title])
        Makie.display(ctx[:figure])
    else
        ctx[:data][]=(g=grid,x=ctx[:xplane],y=ctx[:yplane],z=ctx[:zplane])
        yieldwait()
    end
    ctx[:figure]
end

function plot!(ctx, ::Type{MakieType}, ::Type{Val{3}}, grid , func)
    
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
        ctx[:scene]=Makie.Scene(scale_plot=false)
        ctx[:data]=Makie.Node((g=grid,f=func,x=ctx[:xplane],y=ctx[:yplane],z=ctx[:zplane],l=ctx[:flevel]))
                                  
        if ctx[:outline]


            ctx[:facedata]=Makie.lift(d->extract_visible_bfaces3D(d.g,
                                                                  xyzmax,
                                                                  primepoints=hcat(xyzmin,xyzmax),
                                                                  Tp=Point3f0,
                                                                  Tf=GLTriangleFace),
                                      ctx[:data])
            ctx[:facemeshes]=Makie.lift(d->[make_mesh(d[1][i],d[2][i]) for i=1:nbregions], ctx[:facedata])
            for i=1:nbregions
                Makie.mesh!(ctx[:scene],Makie.lift(d->d[i], ctx[:facemeshes]),
                color=(frgb(Makie,i,nbregions,pastel=true),ctx[:alpha]),
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
        add_scene!(ctx[:ax],ctx[:scene],title=ctx[:title])
        Makie.display(ctx[:figure])
        
        scene_interaction(ctx[:scene],[:z,:y,:x,:l]) do delta,key
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

