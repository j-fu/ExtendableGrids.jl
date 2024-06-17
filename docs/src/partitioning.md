# Grid partitioning


## Access to partitioning data

### API calls
```@docs
num_pcolors
num_partitions
num_partitions_per_color
pcolors
pcolor_partitions
partition_cells
partition_nodes
checkpartitioning
```

### Key types for grid access
```@docs
PColorPartitions 
PartitionCells
PartitionNodes
NodePermutation
```

## Partitioning algorithms
```@docs
AbstractPartitioningAlgorithm
partition
TrivialPartitioning
PlainMetisPartitioning
```



## Internal API
```@docs
ExtendableGrids.trivial_partitioning!
ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PColorPartitions})
ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PartitionCells})
ExtendableGrids.partgraph
ExtendableGrids.induce_node_partitioning!
```
