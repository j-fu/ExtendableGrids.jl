"""
    $(TYPEDEF)

Key type describing colors of partitions. These correspond to 
a coloring of the neigborhood graphs of partitions such that 
operations (e.g. FEM assembly) on partitions of a given color can 
be performed in parallel.

`grid[PColorPartitions]` returns an integer vector describing 
the partition colors ("pcolors") of a grid. 
Let `p=grid[PColorPartitions]`. Then all partitions with numbers
`i ∈ p[c]:p[c+1]-1`  have "color" `c`. 
"""
abstract type PColorPartitions <: AbstractGridIntegerArray1D end


"""
    $(TYPEDEF)

Key type describing the cells of a given partition.

`grid[PartitionCells]` returns an integer vector describing 
the cells of a partition  given by its number.
Let `pc=grid[PartitionCells]`. Then all cells with index
`i ∈ pc[p]:pc[p+1]-1`  belong to partition p.
"""
abstract type PartitionCells <: AbstractGridIntegerArray1D end


"""
    $(TYPEDEF)

Key type describing the bondary faces of a given partition.

`grid[PartitionBFaces]` returns an integer vector describing 
the boundary faces of a partition  given by its number.
Let `pc=grid[PartitionCells]`. Then all cells with index
`i ∈ pc[p]:pc[p+1]-1`  belong to partition p.
"""
abstract type PartitionBFaces <: AbstractGridIntegerArray1D end

"""
    $(TYPEDEF)

Key type describing the nodes of a given partition.

`grid[PartitionNodes]` returns an integer vector describing 
the nodes of a partition  given by its number.
Let `pn=grid[PartitionNodes]`. Then all nodes with index
`i ∈ pn[p]:pn[p+1]-1`  belong to partition p.
"""
abstract type PartitionNodes <: AbstractGridIntegerArray1D end

"""
    $(TYPEDEF)

Key type describing the edges of a given partition.

`grid[PartitionEdges]` returns an integer vector describing 
the edges of a partition  given by its number.
Let `pe=grid[PartitionEdges]`. Then all edges with index
`i ∈ pe[p]:pe[p+1]-1`  belong to partition p.
"""
abstract type PartitionEdges <: AbstractGridIntegerArray1D end


"""
    $(TYPEDEF)

Key type describing the permutation of the nodes of a partitioned grid
with respect  to the unpartitioned origin.

If `pgrid` is the partitioned grid and `grid` is the unpartitioned origin,
then 

`pgrid[Coordinates][:,pgrid[NodePermutation]]==grid[Coordinates]`

"""
abstract type NodePermutation <: AbstractGridIntegerArray1D end

"""
$(SIGNATURES)


(internal)
Create trivial partitioning: the whole grid is partition #1 with just one color.
"""
function trivial_partitioning!(grid::ExtendableGrid{Tc,Ti}) where {Tc,Ti}
    grid[PColorPartitions]=[1,2]
    grid[PartitionCells]=Ti[1,num_cells(grid)+1]
    grid[PartitionBFaces]=Ti[1,num_bfaces(grid)+1]
    grid[PartitionNodes]=Ti[1,num_nodes(grid)+1]
    grid
end


"""
    instantiate(grid::ExtendableGrid, ::Type{PColorPartitions})

If not given otherwise, instantiate partition data with trivial partitioning.
"""
function ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PColorPartitions})
    trivial_partitioning!(grid)
    grid[PColorPartitions]
end

"""
    instantiate(grid::ExtendableGrid, ::Type{PartitionCells})

If not given otherwise, instantiate partition data with trivial partitioning.
"""
function ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PartitionCells})
    trivial_partitioning!(grid)
    grid[PartitionCells]
end

"""
    instantiate(grid::ExtendableGrid, ::Type{PartitionBFaces})

If not given otherwise, instantiate partition data with trivial partitioning.
"""
function ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PartitionBFaces})
    trivial_partitioning!(grid)
    grid[PartitionBFaces]
end

"""
    instantiate(grid::ExtendableGrid, ::Type{PartitionNodes})

If not given otherwise, instantiate partition data with trivial partitioning.
"""
function ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PartitionNodes})
    trivial_partitioning!(grid)
    grid[PartitionNodes]
end

"""
    instantiate(grid::ExtendableGrid, ::Type{PartitionEdges})

If not given otherwise, instantiate partition data with trivial partitioning.
"""
function ExtendableGrids.instantiate(grid::ExtendableGrid{Tc,Ti}, ::Type{PartitionEdges}) where {Tc, Ti}
    trivial_partitioning!(grid)
    grid[PartitionEdges]=Ti[1,num_edges(grid)+1]
    grid[PartitionEdges]
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

Return range of all pcolors.
"""
pcolors(grid)=1:num_pcolors(grid)

"""
    $(SIGNATURES)

Return range of partitions for given pcolor.
"""
function pcolor_partitions(grid,color)
    colpart=grid[PColorPartitions]
    @inbounds colpart[color]:colpart[color+1]-1
end


"""
    $(SIGNATURES)

Return range of cells belonging to a given partition.
"""
function partition_cells(grid, part)
    partcells=grid[PartitionCells]
    @inbounds partcells[part]:partcells[part+1]-1
end

"""
    $(SIGNATURES)

Return range of boundary faces belonging to a given partition.
"""
function partition_bfaces(grid, part)
    partbfaces=grid[PartitionBFaces]
    @inbounds partbfaces[part]:partbfaces[part+1]-1
end


"""
    $(SIGNATURES)

Return range of nodes belonging to a given partition.
"""
function partition_nodes(grid, part)
    partnodes=grid[PartitionNodes]
    partnodes[part]:partnodes[part+1]-1
end

"""
    $(SIGNATURES)

Return range of edges belonging to a given partition.
"""
function partition_edges(grid, part)
    partedges=grid[PartitionEdges]
    partedges[part]:partedges[part+1]-1
end


"""
    $(SIGNATURES)

Return a vector containing the number of partitions for each of
the colors of the grid partitioning.
"""
num_partitions_per_color(grid) = [length(pcolor_partitions(grid, col)) for col in pcolors(grid)]

"""
    $(SIGNATURES)

Return a vector containing the number of cells for each of
the colors of the grid partitioning.
"""
function num_cells_per_color(grid)
    [sum( p->length(partition_cells(grid,p)), pcolor_partitions(grid, col)) for col in pcolors(grid)]
end

"""
    $(SIGNATURES)

Return a vector containing the number of nodes for each of
the partitions of the grid partitioning.
"""
function num_nodes_per_partition(grid)
    @show [length(partition_nodes(grid,ipart)) for ipart=1:num_partitions(grid)]
end

"""
    $(SIGNATURES)

Return a vector containing the number of nodes for each of
the partitions of the grid partitioning.
"""
function num_edges_per_partition(grid)
    @show [length(partition_edges(grid,ipart)) for ipart=1:num_partitions(grid)]
end

"""
    check_partitioning(grid; 
                       verbose=true, 
                       cellpartonly=false)

Check correctness of cell partitioning, necessesary for parallel assembly:
- Check if every node belongs to one of the cell partitions
- Check if no node belongs to two cell partitions of the same color at once



If `cellpartonly==false` check correctness of node partitioning necessary
for parallel sparse matrix multiplication and ILU preconditioning
- Check if no node belongs to two node partitions of the same color at once
- Check if no node is a neighbor of nodes from two node partitions of the same color
"""
function check_partitioning(grid::ExtendableGrid{Tc, Ti}; verbose=true, cellpartonly=false) where {Tc, Ti}
    cn=grid[CellNodes]
    ok=true
    partnodes=Vector{Tc}[unique(vec(cn[:,partition_cells(grid,ipart)])) for ipart=1:num_partitions(grid)]

    if verbose
        @info "Check if every node belongs to one of the cell partitions..."
    end
    if length(intersect(vcat(partnodes...), 1:num_nodes(grid))) !=num_nodes(grid)
        if verbose
            @warn "Not all nodes in one of the partitions"
        end
        ok=false
    end
    
    if verbose
        @info "Check if no node belongs to two cell partitions of the same color at once..."
    end
    for color in pcolors(grid)
        for ipart in pcolor_partitions(grid, color)
            for jpart in pcolor_partitions(grid, color)
                if ipart != jpart
                    is=intersect(partnodes[ipart],partnodes[jpart])
                    if length(is)>0
                        if verbose
                            @warn "Found nodes belonging to cell partitions $ipart,$jpart of color $color"
                        end
                        ok=false
                    end
                end
            end
        end
    end
    
    if cellpartonly
        return ok
    end
    if verbose
        @info "Check if no node belongs to two node partitions of the same color at once..."
    end
    pnodes=grid[PartitionNodes]
    partnodes=Vector{Tc}[collect(pnodes[ipart]:pnodes[ipart+1]-1) for ipart=1:num_partitions(grid)]
    for color in pcolors(grid)
        for ipart in pcolor_partitions(grid, color)
            for jpart in pcolor_partitions(grid, color)
                if ipart != jpart
                    is=intersect(partnodes[ipart],partnodes[jpart])
                    if length(is)>0
                        if verbose
                            @warn "Found nodes belonging to both node partitions $ipart,$jpart of color $color"
                        end
                        ok=false
                    end
                end
            end
        end
    end

    if verbose
        @info "Check if no node is a neighbor of nodes from  two node partitions of the same color..."
    end
    nc = asparse(atranspose(cn))
    rv=SparseArrays.getrowval(nc)
    mark=zeros(Int, num_nodes(grid))
    for color in pcolors(grid)
        mark.=0
        for ipart in pcolor_partitions(grid, color)
            for inode in pnodes[ipart]:pnodes[ipart+1]-1
                for j in nzrange(nc,inode)
                    icell=rv[j]
                    for k=1:size(cn,1)
                        mpart=mark[cn[k,icell]]
                        if mpart!=0 && mpart !=ipart
                            if verbose 
                                @warn "Found node $(cn[k,icell]) neigboring to both partitions $ipart,$mpart of color $color"
                            end
                            ok=false
                        end
                        mark[cn[k,icell]]=ipart
                    end
                end
            end
        end
    end
    if verbose && !ok
        throw(ErrorException("Inconsistency in grid partitioning. Errors in assembly and  matrix-vector multiplication may occur."))
    end
    ok
end

"""
    $(SIGNATURES)

(internal)
Create neigbourhood graph for given partitioning.
"""
function partgraph(cellpartitions,ncellpartitions,cellcelladj)
    gr=Graphs.Graph(ncellpartitions)
    for ic = 1:size(cellcelladj,2)
        ir = cellpartitions[ic]
        for k = cellcelladj.colptr[ic]:cellcelladj.colptr[ic+1]-1
            jr = cellpartitions[cellcelladj.rowval[k]]
            if ir != jr
		Graphs.add_edge!(gr,ir,jr)
            end
        end
    end
    gr
end


"""
    $(TYPEDEF)

Abstract super type for partitioning algorithms
"""
abstract type AbstractPartitioningAlgorithm end



"""
    $(TYPEDEF)

Trivial partitioning: all grid cells belong to single partition number 1.
"""
Base.@kwdef struct TrivialPartitioning <: AbstractPartitioningAlgorithm
end



"""
    $(TYPEDEF)

Subdivide grid into `npart` partitions using `Metis.partition` and color the resulting partition neigborhood graph.
This requires to import Metis.jl in order to trigger the corresponding extension.

Parameters: 

$(TYPEDFIELDS)
"""
Base.@kwdef struct PlainMetisPartitioning <: AbstractPartitioningAlgorithm
    "Number of partitions (default: 20)"
    npart::Int=20

    "Induce node partioning (default: true)"
    partition_nodes::Bool=true

    "Induce edge partioning (default: true)"
    partition_edges::Bool=true

    "Keep node permutation vector (default: true)"
    keep_nodepermutation::Bool=true
end


"""
    $(TYPEDEF)

Subdivide grid into `npart` partitions using `Metis.partition` and calculate cell separators
from this partitioning. The initial partitions  gets pcolor 1, and the separator gets pcolor 2.
This is continued recursively with partitioning of the separator.  

Parameters: 

$(TYPEDFIELDS)
"""
Base.@kwdef struct RecursiveMetisPartitioning <: AbstractPartitioningAlgorithm
    "Number of level 0 partitions (default: 4)"
    npart::Int=4

    "Recursion depth (default: 1)"
    maxdepth::Int=1

    "Separator width  (default: 2)"
    separatorwidth::Int=2

    "Induce node partioning (default: true)"
    partition_nodes::Bool=true

    "Induce edge partioning (default: true)"
    partition_edges::Bool=true

    "Keep node permutation vector (default: true)"
    keep_nodepermutation::Bool=true
end



"""
    $(SIGNATURES)

(Internal utility function)
Create cell permutation such that  all cells belonging to one partition
are contiguous and reorder return grid with reordered cells.
"""
function reorder_cells(grid::ExtendableGrid{Tc,Ti}, cellpartitions,ncellpartitions,colpart) where {Tc,Ti}
    ncells=num_cells(grid)
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
    pgrid[PartitionBFaces]=trivial_partitioning(ncellpartitions,num_bfaces(grid))

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

"""
    $(SIGNATURES)

(internal)
Create a trivial partitioning such that all items
fall in the first of nparts
"""
function trivial_partitioning(npart,nitems)
    part=zeros(Int,npart+1)
    part[1]=1
    part[2:end].=nitems+1
    part
end


"""
    $(SIGNATURES)

(internal)
Induce node partitioning from cell partitioning of `grid`.

Node partitioning should support parallel matrix-vector products with `SparseMatrixCSC`.
The current algorithm assumes that nodes get the partition number from the partition
numbers of the cells having this node in common. If these are differnt, the highest
number is taken.

This algorithm does not always fulfill the  condition that
there is no node which is neigbour of nodes from two different partition with the same color.

This situation is detected and corrected by joining respective critical partitions.
This sacrifies parallel efficiency for correctness.
"""
function induce_node_partitioning!(grid::ExtendableGrid{Tc,Ti},cn,nc; trivial=false, keep_nodepermutation=true) where {Tc, Ti}
    partcells=grid[PartitionCells]
    nnodepartitions=length(partcells)-1
    if trivial
        grid[PartitionNodes]=trivial_partitioning(nnodepartitions,num_nodes(grid))
        if keep_nodepermutation
            grid[NodePermutation]=1:num_nodes(grid)
        end
        return grid
    end

    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    bfacenodes=grid[BFaceNodes]
    lnodes=size(cellnodes,1)

    nodepartitions=zeros(Ti,num_nodes(grid))
    for ipart=1:length(partcells)-1
        for icell in partition_cells(grid,ipart)
            for k=1:lnodes
                nodepartitions[cellnodes[k,icell]]=ipart
            end
        end
    end


    partcolors=zeros(Int,num_partitions(grid))
    for col in pcolors(grid)
        for part in pcolor_partitions(grid,col)
            partcolors[part]=col
        end
    end

    # Correct situation where a node is a neighbor of
    # two different partitions of the same color
    # which would lead to clashes in the matrix-vector product
    nn=cn*nc
    rv=SparseArrays.getrowval(nn)
    showinfo=true

    # Repeat several times until ok.
    while true
        # Detect the situation, record  the corresponding
        # pairs of partitions
        idpart=Pair{Int,Int}[]
        for inode=1:num_nodes(grid)
            ipart=nodepartitions[inode]
            icol=partcolors[ipart]
            for j in nzrange(nn,inode)
                jnode=rv[j]
                for m in nzrange(nn,jnode)
                    mnode=rv[m]
                    mpart=nodepartitions[mnode]
                    mcol=partcolors[mpart]
                    if (mpart != ipart) && (mcol == icol)
                        if ipart > mpart
                            push!(idpart, ipart=>mpart)
                        else
                            push!(idpart, mpart=>ipart)
                        end
                    end
                end
            end
        end
        idpart=unique(idpart)
        # We are done if no such case is left
        if length(idpart)==0
            break
        end
        # Show this only once:
        if showinfo
           @info """Renumber node partitions such that no node is neighbor of
                         two different partitions with same color:\n"""
            showinfo=false
        end
        @info "Renumbering: $idpart\n"

        # Re-assign the respective lower partition numbers for the problem cases.
        # The parttions with the higher numbers will be just empty.
        # This renders the critical nodes into the same partition, so they
        # will be accessed from the same parallel task.
        for inode=1:num_nodes(grid)
            part=nodepartitions[inode]
            for id in idpart
                if part==id[1]
                    nodepartitions[inode]=id[2]
                end
            end
        end
    end

    # Create node permutation such that
    # all nodes belonging to one partition
    # are contiguous
    nodeperm=copy(nodepartitions)
    partctr=zeros(Int,nnodepartitions+1)
    partctr[1]=1
    for i=1:nnodepartitions
	partctr[i+1]=partctr[i]+sum(x->x==i, nodepartitions)
    end
    partnodes=copy(partctr)
    for inode=1:num_nodes(grid)
	part=nodepartitions[inode]
	nodeperm[partctr[part]]=inode
	partctr[part]+=1
    end

    # Permute coordinate array
    nodeperm=invperm(nodeperm)
    xcoord=similar(coord)
    xcoord[:,nodeperm].=coord

  
    # Renumber node indices for cells
    xcellnodes=similar(cellnodes)
    for icell=1:num_cells(grid)
        for k=1:lnodes
            xcellnodes[k,icell]=nodeperm[cellnodes[k,icell]]
        end
    end

    # Renumber node indices for bfaces
    xbfacenodes=similar(bfacenodes)
    for ibface=1:num_bfaces(grid)
        for k=1:lnodes-1
            xbfacenodes[k,ibface]=nodeperm[bfacenodes[k,ibface]]
        end
    end


    grid[PartitionNodes]=partnodes
    grid[CellNodes]=xcellnodes
    grid[BFaceNodes]=xbfacenodes
    grid[Coordinates]=xcoord
    if keep_nodepermutation
        grid[NodePermutation]=nodeperm
    end
    grid
end

function induce_edge_partitioning!(grid::ExtendableGrid{Tc,Ti},cn,nc; trivial=false) where {Tc, Ti}
    partcells=grid[PartitionCells]
    nedgepartitions=length(partcells)-1
    if trivial
        grid[PartitionEdges]=trivial_partitioning(nedgepartitions,num_edges(grid))
        return grid
    end
    coord=grid[Coordinates]
    celledges=grid[CellEdges]
    ledges=size(celledges,1)

    edgepartitions=zeros(Ti,num_edges(grid))
    for ipart=1:length(partcells)-1
        for icell in partition_cells(grid,ipart)
            for k=1:ledges
                edgepartitions[celledges[k,icell]]=ipart
            end
        end
    end


    partcolors=zeros(Int,num_partitions(grid))
    for col in pcolors(grid)
        for part in pcolor_partitions(grid,col)
            partcolors[part]=col
        end
    end

    # Create edge permutation such that
    # all edges belonging to one partition
    # are contiguous
    edgeperm=copy(edgepartitions)
    partctr=zeros(Int,nedgepartitions+1)
    partctr[1]=1
    for i=1:nedgepartitions
	partctr[i+1]=partctr[i]+sum(x->x==i, edgepartitions)
    end
    partedges=copy(partctr)
    for iedge=1:num_edges(grid)
	part=edgepartitions[iedge]
	edgeperm[partctr[part]]=iedge
	partctr[part]+=1
    end

    invedgeperm=invperm(edgeperm)
  
    # Renumber edge indices for cells
    xcelledges=similar(celledges)
    for icell=1:num_cells(grid)
        for k=1:ledges
            xcelledges[k,icell]=invedgeperm[celledges[k,icell]]
        end
    end

    grid[PartitionEdges]=partedges
    grid[CellEdges] = xcelledges
    grid[EdgeCells] = grid[EdgeCells][:,edgeperm]
    grid[EdgeNodes] = grid[EdgeNodes][:,edgeperm]

    # xgrid[CellEdgeSigns] is not changed
    # xgrid[EdgeGeometries] is not changed

    grid
end


"""
    partition(grid, alg::AbstractPartitioningAlgorithm)

Partition grid according to `alg`, such that the neigborhood graph
of partitions is colored in such a way, that all partitions with 
a given color can be worked on in parallel.
"""
function partition(grid::ExtendableGrid, alg::AbstractPartitioningAlgorithm)
    dopartition(grid,alg)
end

"""
    $(SIGNATURES)

(Internal utility function)
Core method for partitioning
"""
function dopartition(grid::ExtendableGrid, alg::AbstractPartitioningAlgorithm)
    if isa(alg,PlainMetisPartitioning)
        error("Import Metis.jl to allow Metis based partitioning")
    else
        error("unexpected Partitioning")
    end
end

function dopartition(grid::ExtendableGrid{Tc,Ti}, ::TrivialPartitioning) where {Tc,Ti}
    pgrid=ExtendableGrid{Tc,Ti}()
    for (k,v) in pairs(grid.components)
        pgrid.components[k]=v
    end
    trivial_partitioning!(pgrid)
end
