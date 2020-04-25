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
Currently, only for 1D grids, inner boundaries can be added.
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
    dim=dim_space(grid)
    bfaceregions=grid[BFaceRegions]
    coord=grid[Coordinates]
    
    function isbface(ix)
        for ibface=1:num_bfaces(grid)
            if bfacenodes[1,ibface]==ix
                return ibface
            end
            return 0
        end
    end
    if dim_space(grid)==1
        Ti=eltype(bfacenodes)
        bfacenodes=ElasticArray{Ti,2}(bfacenodes)
        for inode=1:num_nodes(grid)
            x=coord[1,inode]
            if x>xmaskmin[1] && x<xmaskmax[1]
                ibface=isbface(inode)
                if ibface>0
                    bfaceregions[ibface]=ireg
                else
                    ibface=length(bfaceregions)+1
                    push!(bfaceregions,ireg)
                    append!(bfacenodes,[inode])
                end
            end
        end
        grid[BFaceNodes]=Array{Ti,2}(bfacenodes)
    else
        for ibface=1:num_bfaces(grid)
            in_region=true
            for inode=1:num_targets(bfacenodes,ibface)
                ignode=bfacenodes[inode,ibface]
                for idim=1:dim_space(grid)
                    if coord[idim,ignode]<xmaskmin[idim]
                        in_region=false
                    elseif coord[idim,ignode]>xmaskmax[idim]
                        in_region=false
                    end
                end
            end
            if in_region
                bfaceregions[ibface]=ireg
            end
        end
    end

    grid[NumBFaceRegions]=max(num_bfaceregions(grid),ireg)
    return grid
end
