abstract type PColorPartitions <: AbstractGridIntegerArray1D end
abstract type PartitionCells <: AbstractGridIntegerArray1D end

"""
$(SIGNATURES)

Create default partitioning: the whole grid is partition #1.
"""
function default_partitioning!(grid::ExtendableGrid{Tc,Ti}) where {Tc,Ti}
    grid[PColorPartitions]=ones(Ti,2)
    grid[PartitionCells]=Ti[1,num_cells(grid)]
end

function ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PColorPartitions})
    default_partitioning!(grid)
    grid[PColorPartitions]
end

function ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PartitionCells})
    default_partitioning!(grid)
    grid[PartionCells]
end


"""
$(SIGNATURES)

Return number of partition colors.
"""
num_pcolors(grid)=length(grid[PColorPartitions])-1

"""
$(SIGNATURES)

Return number of partitions.
"""
num_partitions(grid)=length(grid[PartitionCells])-1


"""
    $(SIGNATURES)

Return range of colors.
"""
pcolors(grid)=1:num_pcolors(grid)

"""
    $(SIGNATURES)

Return range of partitions for given color.
"""
function pcolor_partitions(grid,color)
    colpart=grid[PColorPartitions]
    colpart[color]:colpart[color+1]-1
end

"""
    $(SIGNATURES)

Return range of cells for given partition.
"""
function partition_cells(grid, part)
    partcells=grid[PartitionCells]
    partcells[part]:partcells[part+1]-1
end

"""
    $(SIGNATURES)

Partition grid into `npart` partitions using Metis.
Return vector of partition numbers per cell and cell-to-cell
adjacency.
"""
function metis_partition(g::ExtendableGrid, npart)
    nn = num_nodes(g)
    ncells = num_cells(g)
    cn = asparse(g[CellNodes])
    nc = asparse(atranspose(g[CellNodes]))
    cc = nc * cn
    mg = Metis.graph(cc)
    p = Metis.partition(mg, npart)
    p, cc
end

"""
    $(SIGNATURES)

Create neigbouring graph for given partition.
"""
function partgraph(p,np,cc)
    gr=Graphs.Graph(np)
    for ic = 1:cc.n
        ir = p[ic]
        for k = cc.colptr[ic]:cc.colptr[ic+1]-1
            jr = p[cc.rowval[k]]
            if ir != jr
		Graphs.add_edge!(gr,ir,jr)
            end
        end
    end
    gr
end




function partition(grid::ExtendableGrid{Tc,Ti}; npart=20) where {Tc,Ti}
    p,cc=metis_partition(grid,npart)
    up=unique(p)
    np=length(up)
    gr=partgraph(p,np,cc)
    col=Graphs.degree_greedy_color(gr)

    # Reorder partitions such that all partitions with the same color
    # are numbered contiguously
    partperm=similar(up)
    colctr=zeros(Int,col.num_colors+1)
    colctr[1]=1
    for i=1:length(colctr)-1
	colctr[i+1]=colctr[i]+sum(x->x==i, col.colors)
    end
    colpart=copy(colctr)
    for ip=1:length(up)
	color=col.colors[ip]
	partperm[ip]=colctr[color]
	colctr[color]+=1
    end
    
    # Renumber cell partition according to new order
    for i=1:length(p)
        p[i]=partperm[p[i]]
    end

    # Create cell permutation such that
    # all cells belonging to one partition
    # are contiguous
    cellperm=copy(p)
    partctr=zeros(Int,np+1)
    partctr[1]=1
    for i=1:length(partctr)-1
	partctr[i+1]=partctr[i]+sum(x->x==i, p)
    end
    partcells=copy(partctr)
    for ic=1:length(p)
	part=p[ic]
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
