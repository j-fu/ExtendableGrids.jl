import Metis

@testset "partitioning"  begin
    for dim in 1:3
        @info "Test partitioning for dim=$(dim)"
        X=0:0.1:10
        if dim==3
            X=0:0.5:10
        end
        
        XXX=(X for i=1:dim)
        grid1=simplexgrid(XXX...)
        grid2=partition(grid1,TrivialPartitioning())
        grid3=partition(grid1,PlainMetisPartitioning(npart=1))
        grid3=partition(grid1,RecursiveMetisPartitioning(npart=1))
        
        
        for grid in (grid1, grid2, grid3)                
            @test num_pcolors(grid)==1
            @test num_partitions(grid)==1
            @test pcolors(grid) == 1:1
            @test pcolor_partitions(grid,1) == 1:1
            @test partition_cells(grid,1) == 1:num_cells(grid)
            @test num_partitions_per_color(grid) == [1]
            @test check_partitioning(grid)
        end
        
        # METIS has different results depending on the OS...
        for npart in [10,20]
            grid4=partition(grid1,PlainMetisPartitioning(npart=npart))
            @test num_pcolors(grid4) > 1
            @test num_partitions(grid4)==npart
            @test pcolors(grid4) |> length  >0
            @test partition_cells(grid4,1) |> length >0 
            @test num_partitions_per_color(grid4) |> length >0
            @test grid4[Coordinates][:,grid4[NodePermutation]]≈grid1[Coordinates]
            @test try
                check_partitioning(grid4)
            catch err
                if dim==3 && npart==20
                    println("Catched partitioning error: Ok")
                    true
                else
                    false
                end
            end
        end
        for npart in [3] 
            grid4=partition(grid1,RecursiveMetisPartitioning(npart=npart))
            @test num_pcolors(grid4) > 1
            @test num_partitions(grid4)>npart
            @test pcolors(grid4) |> length  >0
            @test partition_cells(grid4,1) |> length >0 
            @test num_partitions_per_color(grid4) |> length >0
            @test grid4[Coordinates][:,grid4[NodePermutation]]≈grid1[Coordinates]
            @test check_partitioning(grid4)
        end
    end
end
