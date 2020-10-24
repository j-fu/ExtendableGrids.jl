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
end


"""
$(TYPEDSIGNATURES)

Edit region numbers of grid  boundary facets  via rectangular mask.
For 1D grids, inner boundaries can be added by this method.
"""
function bfacemask!(grid::ExtendableGrid,
                    maskmin::AbstractArray,
                    maskmax::AbstractArray,
                    ireg::Int;
                    tol=1.0e-10)

    xmaskmin=maskmin.-tol
    xmaskmax=maskmax.+tol


    nbfaces=num_bfaces(grid)
    bfacenodes=grid[BFaceNodes]
    Ti=eltype(bfacenodes)
    dim=dim_space(grid)
    bfaceregions=grid[BFaceRegions]
    new_bfacenodes=ElasticArray{Ti,2}(bfacenodes)
    coord=grid[Coordinates]
    
    function isbface(ix)
        for ibface=1:num_bfaces(grid)
            if bfacenodes[1,ibface]==ix
                return ibface
            end
        end
        return 0
    end

    function isbface(ix,iy)
        for ibface=1:num_bfaces(grid)
            if (bfacenodes[1,ibface] == ix && bfacenodes[2,ibface] == iy) ||
                (bfacenodes[1,ibface] == iy && bfacenodes[2,ibface] == ix)
                return ibface
            end
        end
        return 0
    end

    if dim_space(grid)==1
        for inode=1:num_nodes(grid)
            x=coord[1,inode]
            if x>xmaskmin[1] && x<xmaskmax[1]
                ibface=isbface(inode)
                if ibface>0
                    bfaceregions[ibface]=ireg
                else
                    push!(bfaceregions,ireg)
                    append!(new_bfacenodes,[inode])
                end
            end
        end
    else
        edgenodes=grid[EdgeNodes]
        for iedge=1:size(edgenodes,2)
            in_region=true
            for inode=1:num_targets(edgenodes,iedge)
                ignode=edgenodes[inode,iedge]
                for idim=1:dim_space(grid)
                    if coord[idim,ignode]<xmaskmin[idim]
                        in_region=false
                    elseif coord[idim,ignode]>xmaskmax[idim]
                        in_region=false
                    end
                end
            end
            if in_region
                ibface=isbface(edgenodes[1,iedge],edgenodes[2,iedge])
                if ibface>0
                    bfaceregions[ibface]=ireg
                else
                    push!(bfaceregions,ireg)
                    append!(new_bfacenodes,[edgenodes[1,iedge],edgenodes[2,iedge]])
                end
            end
        end
    end
    btype=grid[BFaceGeometries][1]
    grid[BFaceNodes]=Array{Ti,2}(new_bfacenodes)
    grid[BFaceGeometries]=VectorOfConstants(btype,length(bfaceregions))
    grid[NumBFaceRegions]=max(num_bfaceregions(grid),ireg)
    return grid
end
