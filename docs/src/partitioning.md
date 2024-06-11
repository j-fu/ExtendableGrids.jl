# Grid partitioning


```@docs
PColorPartitions 
PartitionCells
AbstractPartitioningAlgorithm
partition
TrivialPartitioning
PlainMetisPartitioning
num_pcolors
num_partitions
num_partitions_per_color
pcolors
pcolor_partitions
partition_cells
```





## Internal API
```@docs
ExtendableGrids.trivial_partitioning!
ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PColorPartitions})
ExtendableGrids.instantiate(grid::ExtendableGrid, ::Type{PartitionCells})
ExtendableGrids.metis_partition
ExtendableGrids.partgraph
```
