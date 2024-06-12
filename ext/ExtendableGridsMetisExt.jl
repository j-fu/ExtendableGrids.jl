module ExtendableGridsMetisExt

import Metis
import ExtendableGrids: partition, partgraph
using ExtendableGrids, Graphs

"""
    metis_partition(grid, npart)

Partition grid into `npart` partitions using Metis.
Return vector of partition numbers per cell and cell-to-cell
adjacency.
"""
function metis_partition(g::ExtendableGrid, npart)
    nn = num_nodes(g)
    ncells = num_cells(g)
    cn = asparse(g[CellNodes])
    nc = asparse(atranspose(g[CellNodes]))
    cellcelladj = nc * cn
    mg = Metis.graph(cellcelladj)
    cellpartitions = Metis.partition(mg, npart)
    cellpartitions, cellcelladj
end

function partition(grid::ExtendableGrid{Tc,Ti}, alg::PlainMetisPartitioning) where {Tc,Ti}
    cellpartitions,cellcelladj=metis_partition(grid,alg.npart)
    ncells=length(cellpartitions)
    uniquecellpartitions=unique(cellpartitions)
    ncellpartitions=length(uniquecellpartitions)
    gr=partgraph(cellpartitions,ncellpartitions,cellcelladj)
    col=Graphs.degree_greedy_color(gr)

    # Reorder partitions such that all partitions with the same color
    # are numbered contiguously
    partperm=similar(uniquecellpartitions)
    colctr=zeros(Int,col.num_colors+1)
    colctr[1]=1
    for i=1:length(colctr)-1
	colctr[i+1]=colctr[i]+sum(x->x==i, col.colors)
    end
    colpart=copy(colctr)
    for ip=1:ncellpartitions
	color=col.colors[ip]
	partperm[ip]=colctr[color]
	colctr[color]+=1
    end
    
    # Renumber cell partition according to new order
    for i=1:ncells
        cellpartitions[i]=partperm[cellpartitions[i]]
    end

    # Create cell permutation such that
    # all cells belonging to one partition
    # are contiguous
    cellperm=copy(cellpartitions)
    partctr=zeros(Int,ncellpartitions+1)
    partctr[1]=1
    for i=1:ncellpartitions
	partctr[i+1]=partctr[i]+sum(x->x==i, cellpartitions)
    end
    partcells=copy(partctr)
    for ic=1:ncells
	part=cellpartitions[ic]
	cellperm[partctr[part]]=ic
	partctr[part]+=1
    end

    pgrid=ExtendableGrid{Tc,Ti}()
    pgrid[CellNodes]=grid[CellNodes][:,cellperm]
    pgrid[CellRegions]=grid[CellRegions][cellperm]
    pgrid[PColorPartitions]=colpart
    pgrid[PartitionCells]=partcells
    
    for key in [Coordinates,
                CellGeometries,
                BFaceNodes,
                BFaceRegions,
                BFaceGeometries,
                CoordinateSystem]
        pgrid[key]=grid[key]
    end
    pgrid
end

end
