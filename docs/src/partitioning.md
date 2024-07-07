# Grid partitioning

!!! compat
    Grid partitioning is an experimental feature. Breaking changes
    in this realm may occur with minor version updates.

The general idea is that all grids created from ExtendableGrids
can be considered to be partitioned such that the neighborhood graph  of
the partitions is colored so that operations (FEM/FVM assembly, sparse matrix-vector multiplication
with SparseMatrixCSC) on different partitions of the same color can be performed in parallel 
without writing conflicts in a multithreading environment.


The default partitioning is trivial: all cells and nodes belong to one partition,
and the resulting trivial neighborhood graph is colored with one color.

## API calls
```@docs
partition
pcolors
pcolor_partitions
partition_cells
partition_bfaces
partition_nodes
partition_edges
num_pcolors
num_partitions
num_partitions_per_color
num_nodes_per_partition
num_edges_per_partition
num_cells_per_color
check_partitioning
ExtendableGrids.induce_edge_partitioning!
```


## Partitioning algorithms
```@docs
AbstractPartitioningAlgorithm
TrivialPartitioning
PlainMetisPartitioning
RecursiveMetisPartitioning
```

## Key types for grid access
```@docs
PColorPartitions 
PartitionCells
PartitionBFaces
PartitionNodes
PartitionEdges
NodePermutation
```

## Internal API
These functions & methods are neither exported nor public.
```@docs
ExtendableGrids.trivial_partitioning!
ExtendableGrids.trivial_partitioning
ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PColorPartitions})
ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PartitionCells})
ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PartitionBFaces})
ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PartitionNodes})
ExtendableGrids.partgraph
ExtendableGrids.dopartition
ExtendableGrids.reorder_cells
ExtendableGrids.induce_node_partitioning!
```
