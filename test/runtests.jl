ENV["MPLBACKEND"]="agg"

using Test, ExtendableGrids, GridVisualize
import PyPlot

@testset "Basic" begin
    function test_geomspace()
        X1=geomspace(2.0,3.0,0.2,0.2)
        X2=collect(2:0.2:3)
        @assert X1 ≈  X2
        X1=geomspace(2.0,3.0,0.1,0.2)
        X2=geomspace(2.0,3.0,0.2,0.1)
        
        
        X1=geomspace(2.0,3.0,1.0e-5,0.2)
        X2=geomspace(2.0,3.0,0.2,1.0e-5)
        (X1[2]-X1[1])  ≈ (X2[end]-X2[end-1]) &&
            (X2[2]-X2[1])  ≈ (X1[end]-X1[end-1])
    end
    @test test_geomspace()

    function test_prepare_edges()
        ## Compared with pdelib; have more of these
        nodes=Matrix([
            0.0 1;
            1 0;
            0 1;
            1 1]')
        
        cells=Matrix([
            1 2 3;
            2 3 4
        ]')
        cellmat=[1,1]
        bfaces=Matrix([ 1 2;
                        1 3;
                        2 4]')
        bfacemat=[1,1]
        
        grid=simplexgrid(nodes,cells,cellmat,bfaces,bfacemat)
        
        grid[CellEdges]==[3 5; 2 4; 1 3] &&       
        grid[EdgeNodes]==[2 3 3 4 4; 1 1 2 2 3] &&
            grid[EdgeCells]==[1 1 1 2 2; 0 0 2 0 0] 
        
    end
    @test test_prepare_edges()
end

@testset "Read/Write" begin
    function testrw(grid)
        ftmp=tempname()
        write(ftmp,grid,format="sg")
        grid1=simplexgrid(ftmp,format="sg")
        seemingly_equal(grid1,grid)
    end
    X=collect(0:0.05:1)
    @test testrw(simplexgrid(X))
    @test testrw(simplexgrid(X,X))
    @test testrw(simplexgrid(X,X,X))
end

function testgrid(grid,testdata)
    (num_nodes(grid),num_cells(grid), num_bfaces(grid))==testdata
end

examples1d=joinpath(@__DIR__,"..","examples","examples1d.jl")
include(examples1d)
examples2d=joinpath(@__DIR__,"..","examples","examples2d.jl")
include(examples2d)
examples3d=joinpath(@__DIR__,"..","examples","examples3d.jl")
include(examples3d)
plotting=joinpath(@__DIR__,"..","examples","plotting.jl")
include(plotting)


@testset "1D" begin
    @test testgrid(interval_from_vector(),(21,20,2))
    @test testgrid(interval_localref(),   (27,26,2))
    @test testgrid(interval_multiregion(),(21,20,3))
    @test testgrid(interval_subgrid(),(51,50,2))
end

@testset "2D" begin
    @test testgrid(rectangle(),(441,800,80))
    @test testgrid(rectangle_localref(),(729, 1352, 104))
    @test testgrid(rectangle_multiregion(),(441,800,80))
    @test testgrid(rectangle_subgrid(),(360, 600, 120))
end

@testset "3D" begin
    @test testgrid(quadrilateral(),(330,1200,440))
end

if !Sys.isapple()
    @testset "plotting examples" begin
        include("../docs/makeplots.jl")
        picdir=mktempdir()
        
        @test makeplot("interval_from_vector",picdir)
        @test makeplot("interval_localref",picdir)
        @test makeplot("interval_multiregion",picdir)
        @test makeplot("interval_subgrid",picdir)
        @test makeplot("rectangle",picdir)
        @test makeplot("rectangle_localref",picdir)
        @test makeplot("rectangle_multiregion",picdir)
        @test makeplot("rectangle_subgrid",picdir)
        @test makeplot("quadrilateral",picdir)        
    end
end
