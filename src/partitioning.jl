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
$(SIGNATURES)


(internal)
Create trivial partitioning: the whole grid is partition #1 with just one color.
"""
function trivial_partitioning!(grid::ExtendableGrid{Tc,Ti}) where {Tc,Ti}
    grid[PColorPartitions]=[1,2]
    grid[PartitionCells]=Ti[1,num_cells(grid)+1]
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
    instantiate(grid::ExtendableGrid, ::Type{PColorPartitions})

If not given otherwise, instantiate partition data with trivial partitioning.
"""
function ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PartitionCells})
    trivial_partitioning!(grid)
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

Return range of all pcolors.
"""
pcolors(grid)=1:num_pcolors(grid)

"""
    $(SIGNATURES)

Return range of partitions for given pcolor.
"""
function pcolor_partitions(grid,color)
    colpart=grid[PColorPartitions]
    colpart[color]:colpart[color+1]-1
end

"""
    $(SIGNATURES)

Return range of cells belonging to a given partition.
"""
function partition_cells(grid, part)
    partcells=grid[PartitionCells]
    partcells[part]:partcells[part+1]-1
end


"""
    $(SIGNATURES)

Return a vector containing the number of partitions for each of
the colors of the grid partitioning.
"""
function num_partitions_per_color(grid)
    colpart=grid[PColorPartitions]
    nthd=zeros(Int,0)
    for i=1:length(colpart)-1
        push!(nthd,colpart[i+1]-colpart[i])
    end
    nthd
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

function partition(grid::ExtendableGrid{Tc,Ti}, ::TrivialPartitioning) where {Tc,Ti}
    pgrid=ExtendableGrid{Tc,Ti}()
    for (k,v) in pairs(grid.components)
        pgrid.components[k]=v
    end
    trivial_partitioning!(pgrid)
end


"""
    $(TYPEDEF)

Subdivide grid into `npart` partitions using `Metis.partition` and color the resulting partition neigborhood graph.
This requires to import Metis.jl in order to trigger the corresponding extension.

$(TYPEDFIELDS)
"""
Base.@kwdef struct PlainMetisPartitioning <: AbstractPartitioningAlgorithm
    npart::Int=20
end

"""
    $(SIGNATURES)

Partition grid according to `alg`, such that the neigborhood graph
of partitions is colored in such a way, that all partitions with 
a given color can be worked on in parallel.
"""
function partition(grid, alg::AbstractPartitioningAlgorithm)
    if isa(alg,PlainMetisPartitioning)
        error("Import Metis.jl to allow Metis based partitioning")
    else
        error("unexpected Partitioning")
    end
end
