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
        points[i]=Point2f0(coord[1,i],coord[2,i])
    end
    faces=Vector{GLTriangleFace}(undef,nfaces)
    for i=1:nfaces
        faces[i]=TriangleFace(cellnodes[1,i],cellnodes[2,i],cellnodes[3,i])
    end
    Mesh(points,faces)
end



function make_mesh3(grid::ExtendableGrid,zcut0)
    coord=grid[Coordinates]
    zmin=minimum(coord[3,:])
    zmax=maximum(coord[3,:])
    zcut=zmin+zcut0*(zmax-zmin)

    
    cellnodes=grid[CellNodes]
    points=[ Point3f0(coord[:,i]...) for i=1:size(coord,2)]
    faces=Array{NgonFace{3,Int32},1}(undef,0)
    for itet=1:size(cellnodes,2)
        tet=view(cellnodes,:, itet)
        zcoord=[coord[3,tet[i]]-zcut for i=1:4]
        if mapreduce(z->z<=0,*,zcoord)||mapreduce(z->z>=0,*,zcoord)
            continue
        end
        push!(faces,TriangleFace(tet[1],tet[2],tet[3]))
        push!(faces,TriangleFace(tet[1],tet[2],tet[4]))
        push!(faces,TriangleFace(tet[2],tet[3],tet[4]))
        push!(faces,TriangleFace(tet[3],tet[1],tet[4]))
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

function region_bfacesegments(grid::ExtendableGrid,ibreg,zcut0)
    coord=grid[Coordinates]
    zmin=minimum(coord[3,:])
    zmax=maximum(coord[3,:])
    zcut=zmin+zcut0*(zmax-zmin)
    nbfaces=num_bfaces(grid)
    bfacenodes=grid[BFaceNodes]
    bfaceregions=grid[BFaceRegions]
    points=[Point3f0(coord[:,i]...) for i=1:size(coord,2)]
    faces=Array{GeometryBasics.NgonFace{3,Int32},1}(undef, 0)
    for i=1:nbfaces
        if bfaceregions[i]==ibreg
            tri=view(bfacenodes,:, i)
            zcoord=[coord[3,tri[i]]-zcut for i=1:3]
            if !(mapreduce(z->z>=0,*,zcoord))
                push!(faces,TriangleFace(bfacenodes[:,i]...))
            end
        end
    end
    mesh=GeometryBasics.Mesh(points,faces)
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


function plot!(ctx, ::Type{MakieType}, ::Type{Val{3}},grid)
    Makie=ctx[:Plotter]

    nregions=num_cellregions(grid)
    nregions=1

    make_scene!(ctx)

    # ctx[:zslider] = Makie.slider(LinRange(-0.01,1.01,100), raw = true, camera = Makie.campixel!, start = ctx[:zplane])
    # zplane=ctx[:zslider][end][:value]

    meshes=[make_mesh3(grid,ctx[:zplane])]
    nbregions=num_bfaceregions(grid)
    bsegments=[region_bfacesegments(grid,iregion,ctx[:zplane]) for iregion=1:nbregions]
    
    if !haskey(ctx,:meshes)|| ctx[:num_cellregions]!=nregions
        ctx[:num_cellregions]=nregions
        # if ctx[:aspect]>1.0
        #     Makie.scale!(ctx[:scene],ctx[:aspect],1.0)
        # else
        #     Makie.scale!(ctx[:scene],1.0,1.0/ctx[:aspect])
        # end
        ctx[:meshes]=[Makie.Node(meshes[i]) for i=1:nregions]
        ctx[:bsegments]=[Makie.Node(bsegments[i]) for i=1:nbregions]
        trans=ctx[:alpha]<1
        for i=1:nregions
            Makie.mesh!(ctx[:scene],Makie.lift(a->a, ctx[:meshes][1]), color=(frgb(Makie,i,nregions,pastel=true),ctx[:alpha]),transparency=trans)
            if (!trans)
                Makie.wireframe!(ctx[:scene],Makie.lift(a->a, ctx[:meshes][1]),strokecolor=:black)
            end
        end
        for i=1:nbregions
             Makie.mesh!(ctx[:scene],Makie.lift(a->a, ctx[:bsegments][i]) , color=(frgb(Makie,i,nbregions),ctx[:alpha]),transparency=trans)
            if (!trans)
                Makie.wireframe!(ctx[:scene],Makie.lift(a->a, ctx[:bsegments][i]) , strokecolor=:black)
            end
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

