import Metis

@testset "partitioning"  begin
    
    grid1=simplexgrid(0:0.1:10,0:0.1:10)
    grid2=partition(grid1,TrivialPartitioning())
    grid3=partition(grid1,PlainMetisPartitioning(npart=1))
    
    
    for grid in (grid1, grid2, grid3)                
        @test num_pcolors(grid)==1
        @test num_partitions(grid)==1
        @test pcolors(grid) == 1:1
        @test pcolor_partitions(grid,1) == 1:1
        @test partition_cells(grid,1) == 1:num_cells(grid)
        @test num_partitions_per_color(grid) == [1]
        @test checkpartitioning(grid)
    end

    # METIS has different results depending on the OS...
    grid4=partition(grid1,PlainMetisPartitioning(npart=10))
    @test num_pcolors(grid4) > 1
    @test num_partitions(grid4)==10
    @test pcolors(grid4) |> length  >0
    @test partition_cells(grid4,1) |> length >0 
    @test num_partitions_per_color(grid4) |> length >0
    @test grid4[Coordinates][:,grid4[NodePermutation]]≈grid1[Coordinates]
    @test checkpartitioning(grid4)
    
end
