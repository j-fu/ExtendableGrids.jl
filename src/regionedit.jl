"""
$(TYPEDSIGNATURES)

Edit region numbers of grid cells via rectangular mask.
"""
function cellmask!(grid::ExtendableGrid,
                   maskmin::AbstractArray,
                   maskmax::AbstractArray,
                   ireg::Int;
                   tol=1.0e-10)
    xmaskmin=maskmin.-tol
    xmaskmax=maskmax.+tol
    ncells=num_cells(grid)
    cellnodes=grid[CellNodes]
    dim=dim_space(grid)
    cellregions=grid[CellRegions]
    coord=grid[Coordinates]
    for icell=1:ncells
        in_region=true
        for inode=1:num_targets(cellnodes,icell)
            ignode=cellnodes[inode,icell]
            for idim=1:dim
                if coord[idim,ignode]<xmaskmin[idim]
                    in_region=false
                elseif coord[idim,ignode]>xmaskmax[idim]
                    in_region=false
                end
            end
        end
        if in_region
            cellregions[icell]=ireg
        end
    end
    grid[NumCellRegions]=max(num_cellregions(grid),ireg)
    grid
end


"""
$(TYPEDSIGNATURES)

Edit region numbers of grid  boundary facets  via rectangular mask.
If `allow_new` is true (default), new facets are added
"""
function bfacemask!(grid::ExtendableGrid,
                    maskmin::AbstractArray,
                    maskmax::AbstractArray,
                    ireg::Int;
                    allow_new=true,
                    tol=1.0e-10)

    xmaskmin=maskmin.-tol
    xmaskmax=maskmax.+tol


    bfacenodes=grid[BFaceNodes]
    nbfaces=size(bfacenodes,2)



    
    Ti=eltype(bfacenodes)
    dim=dim_space(grid)
    bfaceregions=grid[BFaceRegions]
    coord=grid[Coordinates]

    xfacenodes=bfacenodes
    if allow_new
        new_bfacenodes=ElasticArray{Ti,2}(bfacenodes)
        facenodes=grid[FaceNodes]
        bfacefaces=grid[BFaceFaces]
        nfaces=size(facenodes,2)
        bmark=zeros(Int,nfaces)
        for ibface=1:nbfaces
            bmark[bfacefaces[ibface]]=ibface
        end
        xfacenodes=facenodes
    end
    
    for ixface=1:size(xfacenodes,2)
        in_region=true
        for inode=1:num_targets(xfacenodes,ixface)
            ignode=xfacenodes[inode,ixface]
            for idim=1:dim_space(grid)
                if coord[idim,ignode]<xmaskmin[idim]
                    in_region=false
                elseif coord[idim,ignode]>xmaskmax[idim]
                    in_region=false
                end
            end
        end

        if in_region
            if allow_new
                ibface=bmark[ixface]
                if ibface>0
                    bfaceregions[ibface]=ireg
                else
                    push!(bfaceregions,ireg)
                    @views append!(new_bfacenodes,xfacenodes[:,ixface])
                end
            else
                bfaceregions[ixface]=ireg
            end
        end
    end


    delete!(grid,BFaceFaces)
    btype=grid[BFaceGeometries][1]
    if allow_new
        grid[BFaceNodes]=Array{Ti,2}(new_bfacenodes)
    end
    grid[BFaceGeometries]=VectorOfConstants{ElementGeometries,Int}(btype,length(bfaceregions))
    grid[NumBFaceRegions]=max(num_bfaceregions(grid),ireg)
    return grid
end


"""
    $(TYPEDSIGNATURES)

Edit region numbers of grid  boundary edges via line mask.
This only works for 3D grids.
"""
function bedgemask!(grid::ExtendableGrid,
                    xa::AbstractArray,
                    xb::AbstractArray,
                    ireg::Int;
                    tol=1.0e-10)
    # Masking of boundary edges makes only sense in 3D
    @assert (dim_space(grid) > 2)

    masked = false

    nbedges        = num_bedges(grid)
    bedgenodes     = grid[BEdgeNodes]
    Ti             = eltype(bedgenodes)
    dim            = dim_space(grid)
    bedgeregions   = grid[BEdgeRegions]
    new_bedgenodes = ElasticArray{Ti,2}(bedgenodes)
    coord          = grid[Coordinates]
    Tv             = eltype(coord)

    # length of boundary edge region
    distsq         = sqrt((xa[1]-xb[1])^2 + (xa[2]-xb[2])^2 + (xa[3]-xb[3])^2)

    bedgenodes = grid[BEdgeNodes]
    # loop over boundary edges
    for ibedge = 1:size(bedgenodes, 2)
        in_region = true
        
        #loop over nodes of boundary edge
        for inode=1:num_targets(bedgenodes, ibedge)
            ignode = bedgenodes[inode, ibedge]

            # we compute the distance of the boundary edge node to the endpoints
            # if the sum of the distances is larger (with tolerance) than the length
            # of the boundary region, the point does not lie on the edge
            distxa = sqrt((xa[1]-coord[1,ignode])^2 
                     + (xa[2]-coord[2,ignode])^2 
                     + (xa[3]-coord[3,ignode])^2)
            distxb = sqrt((coord[1,ignode]-xb[1])^2 
                     + (coord[2,ignode]-xb[2])^2 
                     + (coord[3,ignode]-xb[3])^2)
            diff   = distxa + distxb - distsq 
            if (diff > tol)
                in_region = false
                continue
            end
        end
        
        if in_region
            masked = true
            bedgeregions[ibedge] = ireg
        end
    end
    if !masked
        @warn "Couldn't mask any boundary edges for region $(ireg)"
    end

    grid[NumBEdgeRegions]=max(num_bedgeregions(grid),ireg)
    return grid
end
