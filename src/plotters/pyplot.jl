function initialize_plot!(p, ::Type{PyPlotType})
    PyPlot=p.context[:Plotter]
    PyPlot.PyObject(PyPlot.axes3D)# see https://github.com/JuliaPy/PyPlot.jl/issues/351
    if !haskey(p.context,:figure)
        res=p.context[:resolution]
        p.context[:figure]=PyPlot.figure(p.context[:fignumber],figsize=(res[1]/100,res[2]/100),dpi=100)
        for ctx in p.subplots
            ctx[:figure]=p.context[:figure]
        end
    end
    if p.context[:clear]
        p.context[:figure].clf()
    end
end



"""
$(TYPEDSIGNATURES)
Return tridata to be splatted to PyPlot calls
"""
function tridata(grid)
    coord=grid[Coordinates]
    cellnodes=Matrix(grid[CellNodes])
    coord[1,:], coord[2,:],transpose(cellnodes.-1)
end

# Interfaces to Colors/Collorschemes
rgbtuple(c::RGB)=(red(c),green(c),blue(c))
plaincolormap(ctx)=colorschemes[ctx[:colormap]].colors



### 1D grid
function plot!(ctx, ::Type{PyPlotType}, ::Type{Val{1}}, grid)
    PyPlot=ctx[:Plotter]
    ctx[:ax]=ctx[:figure].add_subplot(ctx[:layout]...,ctx[:iplot])
    if ctx[:clear]
        ctx[:ax].cla()
    end
    ax=ctx[:ax]
    fig=ctx[:figure]

    cellregions=grid[CellRegions]
    cellnodes=grid[CellNodes]
    coord=grid[Coordinates]
    ncellregions=grid[NumCellRegions]
    bfacenodes=grid[BFaceNodes]
    bfaceregions=grid[BFaceRegions]
    nbfaceregions=grid[NumBFaceRegions]
    ncellregions=grid[NumCellRegions]
    
    crflag=ones(Bool,ncellregions)
    brflag=ones(Bool,nbfaceregions)
        
    xmin=minimum(coord)
    xmax=maximum(coord)
    h=(xmax-xmin)/40.0
    ax.set_aspect(ctx[:aspect])
    ax.get_yaxis().set_ticks([])
    ax.set_ylim(-5*h,xmax-xmin)
    cmap=region_cmap(ncellregions)

    for icell=1:num_cells(grid)
        ireg=cellregions[icell]
        label = crflag[ireg] ? "cellregion $(ireg)" : ""
        crflag[ireg]=false
        
        x1=coord[1,cellnodes[1,icell]]
        x2=coord[1,cellnodes[2,icell]]
        ax.plot([x1,x1],[-h,h],linewidth=0.5,color="k",label="")
        ax.plot([x2,x2],[-h,h],linewidth=0.5,color="k",label="")
        ax.plot([x1,x2],[0,0],linewidth=3.0,color=rgbtuple(cmap[cellregions[icell]]),label=label)
    end
    
    cmap=bregion_cmap(ncellregions)
    for ibface=1:num_bfaces(grid)
        ireg=bfaceregions[ibface]
        if ireg >0
            label = brflag[ireg] ? "boundary $(ireg)" : ""
            brflag[ireg]=false
            x1=coord[1,bfacenodes[1,ibface]]
            ax.plot([x1,x1],[-2*h,2*h],linewidth=3.0,color=rgbtuple(cmap[ireg]),label=label)
        end
    end
    if ctx[:legend]
        ax.legend()
    end
    ctx[:figure]
end



### 2D grid
function plot!(ctx, ::Type{PyPlotType}, ::Type{Val{2}},grid)
    PyPlot=ctx[:Plotter]
    ctx[:ax]=ctx[:figure].add_subplot(ctx[:layout]...,ctx[:iplot])
    if ctx[:clear]
        ctx[:ax].cla()
    end
    ax=ctx[:ax]
    fig=ctx[:figure]
    cellregions=grid[CellRegions]
    cellnodes=grid[CellNodes]
    coord=grid[Coordinates]
    ncellregions=grid[NumCellRegions]
    nbfaceregions=grid[NumBFaceRegions]
    ncellregions=grid[NumCellRegions]
    if nbfaceregions>0
        bfacenodes=grid[BFaceNodes]
        bfaceregions=grid[BFaceRegions]
    end
    
    crflag=ones(Bool,ncellregions)
    brflag=ones(Bool,nbfaceregions)
    ax.set_aspect(ctx[:aspect])
    tridat=tridata(grid)
    cmap=region_cmap(ncellregions)
    cdata=ax.tripcolor(tridat...,facecolors=grid[CellRegions],cmap=PyPlot.ColorMap(cmap,length(cmap)))
    cbar=fig.colorbar(cdata,ax=ax,ticks=collect(1:ncellregions))
    if ctx[:edges]
        ax.triplot(tridat...,color="k",linewidth=0.5)
    end


    if nbfaceregions>0 
        cmap=bregion_cmap(nbfaceregions)
        # see https://gist.github.com/gizmaa/7214002
        xc=[coord[:,bfacenodes[1,i]] for i=1:num_sources(bfacenodes)]
        yc=[coord[:,bfacenodes[2,i]] for i=1:num_sources(bfacenodes)]
        rgb=[rgbtuple(cmap[bfaceregions[i]]) for i=1:length(bfaceregions)]
        ax.add_collection(PyPlot.matplotlib.collections.LineCollection(collect(zip(xc,yc)),colors=rgb,linewidth=3))
        for i=1:nbfaceregions
            ax.plot(coord[:,1], coord[:,1],label="b_$(i)", color=rgbtuple(cmap[i]))
        end
    end
    if ctx[:legend]
        ax.legend(loc=ctx[:legend_location])
    end
    ctx[:figure]
end


### 3D Grid
function plot!(ctx, ::Type{PyPlotType}, ::Type{Val{3}},grid)
    # See https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html

    PyPlot=ctx[:Plotter]
    ctx[:ax]=ctx[:figure].add_subplot(ctx[:layout]...,ctx[:iplot],projection="3d")
    ax=ctx[:ax]
    fig=ctx[:figure]
    
    nregions=num_cellregions(grid)
    nbregions=num_bfaceregions(grid)

    xyzmin=zeros(3)
    xyzmax=ones(3)
    coord=grid[Coordinates]
    @views for idim=1:3
        xyzmin[idim]=minimum(coord[idim,:])
        xyzmax[idim]=maximum(coord[idim,:])
    end

    ax.set_xlim3d(xyzmin[1],xyzmax[1])
    ax.set_ylim3d(xyzmin[2],xyzmax[2])
    ax.set_zlim3d(xyzmin[3],xyzmax[3])
    ax.view_init(ctx[:elev],ctx[:azim])

    cmap=region_cmap(nregions)
    bcmap=bregion_cmap(nbregions)

    xyzcut=[ctx[:xplane],ctx[:yplane],ctx[:zplane]]
    
    if ctx[:interior]
        regpoints0,regfacets0=extract_visible_cells3D(grid,
                                                      xyzcut,
                                                      primepoints=hcat(xyzmin,xyzmax)
                                                      )
        regfacets=[reshape(reinterpret(Int32,regfacets0[i]),(3,length(regfacets0[i]))) for i=1:nregions]
        regpoints=[reshape(reinterpret(Float32,regpoints0[i]),(3,length(regpoints0[i]))) for i=1:nregions]

        for ireg=1:nregions
            if size(regfacets[ireg],2)>0
                ax.plot_trisurf(regpoints[ireg][1,:],regpoints[ireg][2,:],transpose(regfacets[ireg].-1),regpoints[ireg][3,:],
                                color=rgbtuple(cmap[ireg]),edgecolors=:black,linewidth=0.5)
            end
        end
    end

    
    
    bregpoints0,bregfacets0=extract_visible_bfaces3D(grid,
                                                     xyzcut,
                                                     primepoints=hcat(xyzmin,xyzmax)
                                                     )
    bregfacets=[reshape(reinterpret(Int32,bregfacets0[i]),(3,length(bregfacets0[i]))) for i=1:nbregions]
    bregpoints=[reshape(reinterpret(Float32,bregpoints0[i]),(3,length(bregpoints0[i]))) for i=1:nbregions]
    for ireg=1:nbregions
        if size(bregfacets[ireg],2)>0
            ax.plot_trisurf(bregpoints[ireg][1,:],bregpoints[ireg][2,:],transpose(bregfacets[ireg].-1),bregpoints[ireg][3,:],
                            color=rgbtuple(bcmap[ireg]),edgecolors=:black,linewidth=0.5)
        end
    end
    



    
    if ctx[:legend]
        ax.legend(loc=ctx[:legend_location])
    end
    ctx[:figure]
end



### 1D Function
function plot!(ctx, ::Type{PyPlotType}, ::Type{Val{1}},grid, func)
    PyPlot=ctx[:Plotter]
    ctx[:ax]=ctx[:figure].add_subplot(ctx[:layout]...,ctx[:iplot])
    if ctx[:clear]
        ctx[:ax].cla()
    end
    ax=ctx[:ax]
    fig=ctx[:figure]
    
    cellnodes=grid[CellNodes]
    coord=grid[Coordinates]
    for icell=1:num_cells(grid)
        i1=cellnodes[1,icell]
        i2=cellnodes[2,icell]
        x1=coord[1,i1]
        x2=coord[1,i2]
        if icell==1 && ctx[:label] !=""
            ax.plot([x1,x2],[func[i1],func[i2]],color=ctx[:color],label=ctx[:label])
        else
            ax.plot([x1,x2],[func[i1],func[i2]],color=ctx[:color])
        end                
    end
    if ctx[:legend]
        ax.legend(loc=ctx[:legend_location])
    end
    if ctx[:axisgrid]
        ax.grid()
    end
    ctx[:figure]
end

### 2D Function
function plot!(ctx, ::Type{PyPlotType}, ::Type{Val{2}},grid, func)
    PyPlot=ctx[:Plotter]
    ctx[:ax]=ctx[:figure].add_subplot(ctx[:layout]...,ctx[:iplot])
    if ctx[:clear]
        ctx[:ax].cla()
    end
    
    ax=ctx[:ax]
    fig=ctx[:figure]
    ax.set_aspect(ctx[:aspect])
    umin=minimum(func)
    umax=maximum(func)
    colorlevels=collect(umin:(umax-umin)/(ctx[:colorlevels]-1):umax)
    isolines=collect(umin:(umax-umin)/(ctx[:isolines]-1):umax)
    if !haskey(ctx,:grid) || !seemingly_equal(ctx[:grid],grid)
        ctx[:grid]=grid
        ctx[:tridata]=tridata(grid)
    end
    cnt=ax.tricontourf(ctx[:tridata]...,func;levels=colorlevels,cmap=PyPlot.ColorMap(plaincolormap(ctx)))
    for c in cnt.collections
        c.set_edgecolor("face")
    end
    if ctx[:colorbar]
        fig.colorbar(cnt,ax=ax,ticks=isolines,boundaries=colorlevels)
    end
    ax.tricontour(ctx[:tridata]...,func,colors="k",levels=isolines)
    ctx[:figure]
end

function plot!(ctx, T::Type{PyPlotType}, ::Type{Val{3}},grid,func)

    PyPlot=ctx[:Plotter]
    ctx[:ax]=ctx[:figure].add_subplot(ctx[:layout]...,ctx[:iplot],projection="3d")
    ax=ctx[:ax]
    fig=ctx[:figure]
    
    nregions=num_cellregions(grid)
    nbregions=num_bfaceregions(grid)

    xyzmin=zeros(3)
    xyzmax=ones(3)
    coord=grid[Coordinates]
    @views for idim=1:3
        xyzmin[idim]=minimum(coord[idim,:])
        xyzmax[idim]=maximum(coord[idim,:])
    end
    xyzcut=[ctx[:xplane],ctx[:yplane],ctx[:zplane]]
    fminmax=extrema(func)
    
    ctx[:xplane]=max(xyzmin[1],min(xyzmax[1],ctx[:xplane]) )
    ctx[:yplane]=max(xyzmin[2],min(xyzmax[2],ctx[:yplane]) )
    ctx[:zplane]=max(xyzmin[3],min(xyzmax[3],ctx[:zplane]) )
    ctx[:flevel]=max(fminmax[1],min(fminmax[2],ctx[:flevel]))

    makeplanes(x,y,z)=[[1,0,0,-x], 
                       [0,1,0,-y], 
                       [0,0,1,-z]]

    ccoord0,faces0,values=marching_tetrahedra(grid,func,makeplanes(ctx[:xplane],ctx[:yplane],ctx[:zplane]),[ctx[:flevel]])

    faces=reshape(reinterpret(Int32,faces0),(3,length(faces0)))
    ccoord=reshape(reinterpret(Float32,ccoord0),(3,length(ccoord0)))

    nfaces=size(faces,2)
    if nfaces>0
        colors=zeros(nfaces)
        for i=1:nfaces
            colors[i]=(values[faces[1,i]]+values[faces[2,i]]+values[faces[3,i]])/3
        end
        # thx, https://stackoverflow.com/a/24229480/8922290 
        collec=ctx[:ax].plot_trisurf(ccoord[1,:],ccoord[2,:],transpose(faces.-1),ccoord[3,:],
                                     cmap=PyPlot.ColorMap(plaincolormap(ctx)))
        collec.set_array(colors)
        collec.autoscale()
    end

    ax.set_xlim3d(xyzmin[1],xyzmax[1])
    ax.set_ylim3d(xyzmin[2],xyzmax[2])
    ax.set_zlim3d(xyzmin[3],xyzmax[3])
    ax.view_init(ctx[:elev],ctx[:azim])
    
    
    if ctx[:legend]
        ax.legend(loc=ctx[:legend_location])
    end
    ctx[:figure]
end

