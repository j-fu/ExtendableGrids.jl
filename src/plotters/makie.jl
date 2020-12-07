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

function make_mesh(grid::ExtendableGrid)
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    npoints=num_nodes(grid)
    nfaces=num_cells(grid)
    
    points=Vector{Point2f0}(undef,npoints)
    for i=1:npoints
        points[i]=Point2f0(coord[1,i],coobrd[2,i])
    end
    faces=Vector{GLTriangleFace}(undef,nfaces)
    for i=1:nfaces
        faces[i]=TriangleFace(cellnodes[1,i],cellnodes[2,i],cellnodes[3,i])
    end
    Mesh(points,faces)
end




function make_mesh(grid::ExtendableGrid, elevation; elevation_factor=1.0)
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    npoints=num_nodes(grid)
    nfaces=num_cells(grid)

    points=Vector{Point3f0}(undef,npoints)
    for i=1:npoints
        points[i]=Point3f0(coord[1,i],coord[2,i],elevation[i]*elevation_factor)
    end
    faces=Vector{GLTriangleFace}(undef,nfaces)
    for i=1:nfaces
        faces[i]=TriangleFace(cellnodes[1,i],cellnodes[2,i],cellnodes[3,i])
    end
    Mesh(points,faces)
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



function cutvalues(coord,xyz)
    xyzcut=[0.0,0.0,0.0]
    @views for idim=1:3
        min=minimum(coord[idim,:])
        max=maximum(coord[idim,:])
        xyzcut[idim]=min+xyz[idim]*(max-min)
    end
    xyzcut
end

function take(coord,simplex,xyzcut)
    nnodes=length(simplex)
    
    # for idim=1:3
    #     cutcoord=[coord[idim,simplex[inode]]-xyzcut[idim] for inode=1:nnodes]
    #     result[idim]=!(mapreduce(a->a<=0,*,cutcoord)|| mapreduce(a->a>=0,*,cutcoord))
    # end
    cutcoord=[]
    for idim=1:3
        for inode=1:nnodes
            push!(cutcoord,coord[idim,simplex[inode]]-xyzcut[idim])
        end
    end
    mapreduce(a->a<=0,*,cutcoord)|| mapreduce(a->a>=0,*,cutcoord)
end

function btake(coord,simplex,xyzcut)
    nnodes=length(simplex)
    for idim=1:3
        cutcoord=[coord[idim,simplex[inode]]-xyzcut[idim] for inode=1:nnodes]
        if !mapreduce(a->a<=0,*,cutcoord)
            return false
        end
    end
    return true
end


function region_bfacesegments(grid::ExtendableGrid,ibreg,xyz)
    coord=grid[Coordinates]
    xyzcut=cutvalues(coord,xyz)
    nbfaces=num_bfaces(grid)
    bfacenodes=grid[BFaceNodes]
    bfaceregions=grid[BFaceRegions]
    points=[Point3f0(coord[:,i]...) for i=1:size(coord,2)]
    faces=Array{GeometryBasics.NgonFace{3,Int32},1}(undef, 0)
    for i=1:nbfaces
        if bfaceregions[i]==ibreg
            tri=view(bfacenodes,:, i)
            if btake(coord,tri,xyzcut)
                push!(faces,TriangleFace(bfacenodes[:,i]...))
            end
        end
    end
    mesh=GeometryBasics.Mesh(points,faces)
end


function make_mesh3(grid::ExtendableGrid,iregion,xyz)
    coord=grid[Coordinates]
    xyzcut=cutvalues(coord,xyz)
    
    cellnodes=grid[CellNodes]
    cellregions=grid[CellRegions]
    points=[ Point3f0(coord[:,i]...) for i=1:size(coord,2)]
    faces=Array{NgonFace{3,Int32},1}(undef,0)
    for itet=1:size(cellnodes,2)
        if cellregions[itet]==iregion
            tet=view(cellnodes,:, itet)
            if take(coord,tet,xyzcut)
                push!(faces,TriangleFace(tet[1],tet[2],tet[3]))
                push!(faces,TriangleFace(tet[1],tet[2],tet[4]))
                push!(faces,TriangleFace(tet[2],tet[3],tet[4]))
                push!(faces,TriangleFace(tet[3],tet[1],tet[4]))
            end
        end
    end
    
    Mesh(points,faces)
end

function plot!(ctx, ::Type{MakieType}, ::Type{Val{3}},grid)
    Makie=ctx[:Plotter]
    if !haskey(ctx,:fullscene)
        ctx[:xslider] = Makie.slider(LinRange(-0.02,1.02,105),start = ctx[:xplane],
                                     camera=Makie.campixel!,
                                     buttoncolor=:gray,
                                     valueprinter= x->@sprintf("x=%.0f%%",100*x),
                                     raw=false)
        ctx[:yslider] = Makie.slider(LinRange(-0.02,1.02,105),start = ctx[:yplane],
                                     camera=Makie.campixel!,
                                     buttoncolor=:gray,
                                     valueprinter= x->@sprintf("y=%.0f%%",100*x),
                                     raw=false)
        ctx[:zslider] = Makie.slider(LinRange(-0.02,1.02,105),start = ctx[:zplane],
                                     camera=Makie.campixel!,
                                     buttoncolor=:gray,
                                     valueprinter= x->@sprintf("z=%.0f%%",100*x),
                                     raw=false)
    end

    

    function makeplot(x,y,z;renew=false)
        alpha=ctx[:alpha]
        trans=alpha<1
        nregions=num_cellregions(grid)
        nbregions=num_bfaceregions(grid)
        
        if (!haskey(ctx,:scene)||
            ctx[:num_cellregions]!=nregions ||
            ctx[:num_bfaceregions]!=nbregions)

            ctx[:scene] = Makie.Scene(scale_plot=false)

            # if ctx[:aspect]>1.0
            #     Makie.scale!(ctx[:scene],ctx[:aspect],1.0)
            # else
            #     Makie.scale!(ctx[:scene],1.0,1.0/ctx[:aspect])
            # end

            ctx[:num_cellregions]=nregions
            ctx[:num_bfaceregions]=nbregions
            ctx[:meshes]=[Makie.Node(make_mesh3(grid,i,[x,y,z])) for i=1:nregions]
            ctx[:bsegments]=[Makie.Node(region_bfacesegments(grid,i,[x,y,z])) for i=1:nbregions]

            if ctx[:interior]
                for i=1:nregions
                    Makie.mesh!(ctx[:scene],Makie.lift(a->a, ctx[:meshes][i]), color=(frgb(Makie,i,nregions,pastel=true),alpha),transparency=trans)
                    if (!trans)
                        Makie.wireframe!(ctx[:scene],Makie.lift(a->a, ctx[:meshes][i]),strokecolor=:black)
                    end
                end
            end
            
            for i=1:nbregions
                Makie.mesh!(ctx[:scene],Makie.lift(a->a, ctx[:bsegments][i]) , color=(frgb(Makie,i,nbregions),alpha),transparency=trans)
                if (!trans)
                    Makie.wireframe!(ctx[:scene],Makie.lift(a->a, ctx[:bsegments][i]) , strokecolor=:black)
                end
            end
            ctx[:fullscene]=Makie.vbox(ctx[:scene],Makie.hbox(ctx[:zslider],ctx[:yslider],ctx[:xslider]))
            Makie.display(ctx[:fullscene])
        else
            if ctx[:interior]
                for i=1:nregions
                    ctx[:meshes][i][]=make_mesh3(grid,i,[x,y,z])
                end
            end
            for i=1:nbregions
                ctx[:bsegments][i][]=region_bfacesegments(grid,i,[x,y,z])
            end
        end
    end
    
    xplane=ctx[:xslider][end][:value]
    yplane=ctx[:yslider][end][:value]
    zplane=ctx[:zslider][end][:value]
    makeplot(xplane[],yplane[],zplane[])
    onany(makeplot,xplane,yplane,zplane)


    # on(ctx[:fullscene].events.unicode_input) do button
    #     if button==['i']
    #         ctx[:interior]=!ctx[:interior]
    #         @show ctx[:interior]
    #         makeplot(xplane[],yplane[],zplane[],renew=true)
    #     end
    # end
            
    
    if ctx[:show]
        Makie.update!(ctx[:fullscene])
    end
    
    ctx[:fullscene]
end



function plot!(ctx, ::Type{MakieType}, ::Type{Val{2}},grid, func)
    Makie=ctx[:Plotter]
    make_scene!(ctx)
    if ctx[:elevation]
        mesh=make_mesh(grid,func,elevation_factor=ctx[:elevation_factor])
    else
        mesh=make_mesh(grid)
    end        
    if !haskey(ctx,:meshnode)
        # if ctx[:aspect]>1.0
        #     Makie.scale!(ctx[:scene],ctx[:aspect],1.0)
        # else
        #     Makie.scale!(ctx[:scene],1.0,1.0/ctx[:aspect])
        # end
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

