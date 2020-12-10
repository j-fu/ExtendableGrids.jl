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


function make_sliders(ctx,xyzmin,xyzmax)
    Makie=ctx[:Plotter]
    scaled_range(idim)=LinRange(-0.02,1.02,105)*(xyzmax[idim]-xyzmin[idim]).+xyzmin[idim]
    ctx[:xslider] = Makie.slider(scaled_range(1),start=ctx[:xplane],
                                 camera=Makie.campixel!,
                                 buttoncolor=:gray,
                                 valueprinter= x->@sprintf("x=%.2g",x),
                                 raw=false)
    ctx[:yslider] = Makie.slider(scaled_range(2),start = ctx[:yplane],
                                 camera=Makie.campixel!,
                                 buttoncolor=:gray,
                                 valueprinter= y->@sprintf("y=%.2g",y),
                                 raw=false)
    ctx[:zslider] = Makie.slider(scaled_range(3),start = ctx[:zplane],
                                 camera=Makie.campixel!,
                                 buttoncolor=:gray,
                                 valueprinter= z->@sprintf("z=%.2g",z),
                                 raw=false)
end


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


function plot!(ctx, ::Type{MakieType}, ::Type{Val{3}},grid)
    Makie=ctx[:Plotter]
    coord=grid[Coordinates]
    xyzmin,xyzmax=xyzminmax(grid)
    
    if !haskey(ctx,:fullscene)
        make_sliders(ctx,xyzmin,xyzmax)
    end


    function makeplot(x,y,z)
        alpha=ctx[:alpha]
        trans=alpha<1
        nregions=num_cellregions(grid)
        nbregions=num_bfaceregions(grid)
        coord=grid[Coordinates]
        xyzcut=[x,y,z]
        
        function makemesh(pts,fcs)
            pts=vec(collect(reinterpret(Point3f0,pts)))
            push!(pts,Point3f0(xyzmin...))
            push!(pts,Point3f0(xyzmax...))
            fcs=vec(collect(reinterpret(NgonFace{3,Int32},fcs)))
            Mesh(pts,fcs)
        end

        if (!haskey(ctx,:scene)||
            ctx[:num_cellregions]!=nregions ||
            ctx[:num_bfaceregions]!=nbregions)

            ctx[:scene] = Makie.Scene(scale_plot=false)

            # TODO: allow aspect scaling
            # if ctx[:aspect]>1.0
            #     Makie.scale!(ctx[:scene],ctx[:aspect],1.0)
            # else
            #     Makie.scale!(ctx[:scene],1.0,1.0/ctx[:aspect])
            # end

            
            regpoints,regfacets=extract_visible_cells3D(grid,xyzcut)
            bregpoints,bregfacets=extract_visible_bfaces3D(grid,xyzcut)
            ctx[:num_cellregions]=nregions
            ctx[:num_bfaceregions]=nbregions
            ctx[:meshes]=   [Makie.Node(makemesh(regpoints[iregion],regfacets[iregion])) for iregion=1:nregions]
            ctx[:bsegments]=[Makie.Node(makemesh(bregpoints[ibregion],bregfacets[ibregion])) for ibregion=1:nbregions]
            
            # TODO: use distinguishable colors
            # http://juliagraphics.github.io/Colors.jl/stable/colormapsandcolorscales/#Generating-distinguishable-colors-1
            if ctx[:interior]
                for i=1:nregions
                    Makie.mesh!(ctx[:scene],Makie.lift(a->a, ctx[:meshes][i]), color=(frgb(Makie,i,nregions+nbregions,pastel=true),alpha),transparency=trans)
                    if (!trans)
                        Makie.wireframe!(ctx[:scene],Makie.lift(a->a, ctx[:meshes][i]),strokecolor=:black)
                    end
                end
            end
            for i=1:nbregions
                Makie.mesh!(ctx[:scene],Makie.lift(a->a, ctx[:bsegments][i]) , color=(frgb(Makie,nregions+i,nregions+nbregions),alpha),transparency=trans)
                if (!trans)
                    Makie.wireframe!(ctx[:scene],Makie.lift(a->a, ctx[:bsegments][i]) , strokecolor=:black)
                end
            end


            # TODO: a priori angles aka pyplot3D
            # rect = ctx[:scene]
            # azim=ctx[:azim]
            # elev=ctx[:elev]
            # arr = normalize([cosd(azim/2), 0, sind(azim/2), -sind(azim/2)])
            # Makie.rotate!(rect, Makie.Quaternionf0(arr...))
            
            ctx[:fullscene]=Makie.vbox(ctx[:scene],Makie.hbox(ctx[:zslider],ctx[:yslider],ctx[:xslider]))
            Makie.display(ctx[:fullscene])

        else
            if ctx[:interior]
                regpoints,regfacets=extract_visible_cells3D(grid,xyzcut)
                for i=1:nregions
                    ctx[:meshes][i][]=makemesh(regpoints[i],regfacets[i])
                end
            end
            bregpoints,bregfacets=extract_visible_bfaces3D(grid,xyzcut)
            for i=1:nbregions
                 ctx[:bsegments][i][]=makemesh(bregpoints[i],bregfacets[i])
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





function plot!(ctx, ::Type{MakieType}, ::Type{Val{3}},grid , func)
    Makie=ctx[:Plotter]

    xyzmin,xyzmax=xyzminmax(grid)
    fminmax = extrema(func)
    cmap = Makie.to_colormap(ctx[:colormap])
 

    if !haskey(ctx,:fullscene)
        make_sliders(ctx,xyzmin,xyzmax)
    end
    
    ctx[:grid]=grid
    ctx[:func]=func
    
    function makeplot(x,y,z)
        alpha=ctx[:alpha]
        trans=alpha<1
        nregions=num_cellregions(grid)
        coord=grid[Coordinates]
        xyzcut=[x,y,z]
        function makemesh(pts,fcs)
            pts=vec(collect(reinterpret(Point3f0,pts)))
            push!(pts,Point3f0(xyzmin...))
            push!(pts,Point3f0(xyzmax...))
            fcs=vec(collect(reinterpret(NgonFace{3,Int32},fcs)))
            Mesh(pts,fcs)
        end

        function makemesh(pts,fcs,vals)
            pts=vec(collect(reinterpret(Point3f0,pts)))
            push!(pts,Point3f0(xyzmin...))
            push!(pts,Point3f0(xyzmax...))
            fcs=vec(collect(reinterpret(NgonFace{3,Int32},fcs)))
            push!(vals,0)
            push!(vals,0)
            colors = Makie.AbstractPlotting.interpolated_getindex.((cmap,), vals, (fminmax,))
            GeometryBasics.Mesh(meta(pts, color=colors), fcs)
        end

        if (!haskey(ctx,:scene)||
            ctx[:num_cellregions]!=nregions)

            ctx[:scene] = Makie.Scene(scale_plot=false)

            # TODO: allow aspect scaling
            # if ctx[:aspect]>1.0
            #     Makie.scale!(ctx[:scene],ctx[:aspect],1.0)
            # else
            #     Makie.scale!(ctx[:scene],1.0,1.0/ctx[:aspect])
            # end

            
            regpoints,regfacets=extract_visible_cells3D(grid,xyzcut)
            bpoints,bfacets,bvalues=extract_visible_bfaces3D(grid,func,xyzcut)
            ctx[:num_cellregions]=nregions
            ctx[:meshes]=  [Makie.Node(makemesh(regpoints[iregion],regfacets[iregion])) for iregion=1:nregions]
            ctx[:bmesh]= Makie.Node(makemesh(bpoints,bfacets,bvalues))
            
            
            # TODO: use distinguishable colors
            # http://juliagraphics.github.io/Colors.jl/stable/colormapsandcolorscales/#Generating-distinguishable-colors-1
            if ctx[:interior]
                for i=1:nregions
                    Makie.mesh!(ctx[:scene],Makie.lift(a->a, ctx[:meshes][i]), color=(frgb(Makie,i,nregions,pastel=true),alpha),transparency=trans)
                    if (!trans)
                        Makie.wireframe!(ctx[:scene],Makie.lift(a->a, ctx[:meshes][i]),strokecolor=:black)
                    end
                end
            end
            
            Makie.mesh!(ctx[:scene], Makie.lift(a->a, ctx[:bmesh]))

            # TODO: a priori angles aka pyplot3D
            # rect = ctx[:scene]
            # azim=ctx[:azim]
            # elev=ctx[:elev]
            # arr = normalize([cosd(azim/2), 0, sind(azim/2), -sind(azim/2)])
            # Makie.rotate!(rect, Makie.Quaternionf0(arr...))
            
            ctx[:fullscene]=Makie.vbox(ctx[:scene],Makie.hbox(ctx[:zslider],ctx[:yslider],ctx[:xslider]))
            Makie.display(ctx[:fullscene])

        else
            if ctx[:interior]
                regpoints,regfacets=extract_visible_cells3D(grid,xyzcut)
                for i=1:nregions
                    ctx[:meshes][i][]=makemesh(regpoints[i],regfacets[i])
                end
            end
            bpoints,bfacets,bvalues=extract_visible_bfaces3D(ctx[:grid],ctx[:func],xyzcut)
            ctx[:bmesh][]= makemesh(bpoints,bfacets,bvalues)
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
