# Grid partitioning

!!! compat
    Grid partitioning is an experimental feature. Breaking changes
    in this realm may occur with minor version updates.

## API calls
```@docs
partition
num_pcolors
num_partitions
num_partitions_per_color
num_cells_per_color
pcolors
pcolor_partitions
partition_cells
partition_nodes
check_partitioning
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
PartitionNodes
NodePermutation
```

## Internal API
```@docs
ExtendableGrids.trivial_partitioning!
ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PColorPartitions})
ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PartitionCells})
ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PartitionNodes})
ExtendableGrids.partgraph
ExtendableGrids.reorder_cells
ExtendableGrids.induce_node_partitioning!
```
