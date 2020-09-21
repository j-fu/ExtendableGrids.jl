function initialize_context!(ctx::PlotterContext,::Type{PlotsType})
    ctx
end
function prepare_plot!(ctx)
    Plots=ctx[:Plotter]
    if !haskey(ctx,:plot) || ctx[:clear]
        ctx[:plot]=Plots.plot(size=ctx[:resolution])
    end
    ctx
end

function plot!(ctx, ::Type{PlotsType}, ::Type{Val{1}},grid)
    Plots=ctx[:Plotter]
    prepare_plot!(ctx)
    
    p=ctx[:plot]
    
    cellregions=grid[CellRegions]
    cellnodes=grid[CellNodes]
    coord=grid[Coordinates]
    ncellregions=grid[NumCellRegions]
    bfacenodes=grid[BFaceNodes]
    bfaceregions=grid[BFaceRegions]
    nbfaceregions=grid[NumBFaceRegions]

    xmin=minimum(coord)
    xmax=maximum(coord)
    h=(xmax-xmin)/20.0
    
    for icell=1:num_cells(grid)
        rgb=frgb(Plots,cellregions[icell],ncellregions)
        x1=coord[1,cellnodes[1,icell]]
        x2=coord[1,cellnodes[2,icell]]
        Plots.plot!(p,[x1,x1],[-h,h],linewidth=0.5,color=:black,label="")
        Plots.plot!(p,[x2,x2],[-h,h],linewidth=0.5,color=:black,label="")
        Plots.plot!(p,[x1,x2],[0,0],linewidth=3.0,color=rgb,label="")
    end
    
    for ibface=1:num_bfaces(grid)
        if bfaceregions[ibface]>0
            rgb=frgb(Plots,bfaceregions[ibface],nbfaceregions)
            x1=coord[1,bfacenodes[1,ibface]]
            Plots.plot!(p,[x1,x1],[-2*h,2*h],linewidth=3.0,color=rgb,label="")
        end
    end
    
    if ctx[:show]
        Plots.gui(p)
    end
    p
end

function plot!(ctx, ::Type{PlotsType}, ::Type{Val{2}},grid)
    Plots=ctx[:Plotter]
    prepare_plot!(ctx)
    p=ctx[:plot]
    
    cellregions=grid[CellRegions]
    cellnodes=grid[CellNodes]
    coord=grid[Coordinates]
    ncellregions=grid[NumCellRegions]
    bfacenodes=grid[BFaceNodes]
    bfaceregions=grid[BFaceRegions]
    nbfaceregions=grid[NumBFaceRegions]

    for icell=1:num_cells(grid)
        rgb=frgb(Plots,cellregions[icell],ncellregions,pastel=true)
        inode1=cellnodes[1,icell]
        inode2=cellnodes[2,icell]
        inode3=cellnodes[3,icell]
        # https://github.com/JuliaPlots/Plots.jl/issues/605
        tri=Plots.Shape([coord[1,inode1],coord[1,inode2], coord[1,inode3]],[coord[2,inode1],coord[2,inode2],coord[2,inode3]])
        Plots.plot!(p,tri,color=rgb,label="")
    end
    for icell=1:num_cells(grid)
        inode1=cellnodes[1,icell]
        inode2=cellnodes[2,icell]
        inode3=cellnodes[3,icell]
        Plots.plot!(p, [coord[1,inode1],coord[1,inode2]],[coord[2,inode1],coord[2,inode2]]  ,linewidth=0.5,color=:black,label="")
        Plots.plot!(p, [coord[1,inode1],coord[1,inode3]],[coord[2,inode1],coord[2,inode3]]  ,linewidth=0.5,color=:black,label="")
        Plots.plot!(p, [coord[1,inode2],coord[1,inode3]],[coord[2,inode2],coord[2,inode3]]  ,linewidth=0.5,color=:black,label="")
    end
    for ibface=1:num_bfaces(grid)
        rgb=frgb(Plots,bfaceregions[ibface],nbfaceregions)
        inode1=bfacenodes[1,ibface]
        inode2=bfacenodes[2,ibface]
        Plots.plot!(p,[coord[1,inode1],coord[1,inode2]],[coord[2,inode1],coord[2,inode2]]  ,linewidth=5,color=rgb,label="")
    end
    if ctx[:show]
        Plots.gui(p)
    end
    p
end

function plot!(ctx, ::Type{PlotsType}, ::Type{Val{1}},grid, func)
    Plots=ctx[:Plotter]
    prepare_plot!(ctx)
    p=ctx[:plot]
    cellnodes=grid[CellNodes]
    coord=grid[Coordinates]
    color=ctx[:color]
    for icell=1:num_cells(grid)
        i1=cellnodes[1,icell]
        i2=cellnodes[2,icell]
        x1=coord[1,i1]
        x2=coord[1,i2]
        if icell==1 && ctx[:label] !=""
            Plots.plot!(p,[x1,x2],[func[i1],func[i2]],linecolor=Plots.RGB(color...),label=ctx[:label])
        else
            Plots.plot!(p,[x1,x2],[func[i1],func[i2]],linecolor=Plots.RGB(color...),label="")
        end                
    end
    if ctx[:show]
        Plots.gui(p)
    end
    p
end



"""
$(TYPEDSIGNATURES)
Return rectangular grid data + function to be splatted into Plots calls
"""
function rectdata(grid,U)
    if dim_grid(grid)==1 && haskey(grid,XCoordinates) 
        return grid[XCoordinates],U
    end
    if dim_grid(grid)==2 && haskey(grid,XCoordinates) && haskey(grid,YCoordinates)
        X=grid[XCoordinates]
        Y=grid[YCoordinates]
        return X,Y,transpose(reshape(U,length(X),length(Y)))
    end
    error("no rectdata on grid")
    nothing
end


function plot!(ctx, ::Type{PlotsType}, ::Type{Val{2}},grid, func)
    Plots=ctx[:Plotter]
    prepare_plot!(ctx)
    p=ctx[:plot]

    rdata=rectdata(grid,func)
    if rdata==nothing
        @error "2D tricontour not available for Plots, see e.g. https://github.com/JuliaPlots/Plots.jl/issues/392"
        return p
    end
    Plots.contourf!(p,rdata...,aspect_ratio=:equal)
    if ctx[:show]
        Plots.gui(p)
    end
    p
end
