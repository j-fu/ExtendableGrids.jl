ENV["MPLBACKEND"]="agg"

using Test, ExtendableGrids, GridVisualize
import PyPlot



@testset "Geomspace" begin
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

    function test_geomspace1(n,h0,h1)
        for i=1:n
            X=geomspace(0,1,h0*rand(),h1*rand())
        end
        true
    end
    
    @test test_geomspace()
    @test test_geomspace1(100,0.01,0.01)
    @test test_geomspace1(100,0.01,0.1)
    @test test_geomspace1(100,0.1,0.01)
    @test test_geomspace1(100,0.001,0.1)


    @test length(geomspace(0,1,0.01,0.1))==27
    @test length(geomspace(0,1,0.001,0.1))==47
    @test length(geomspace(0,1,0.0001,0.1))==68
    @test length(geomspace(0,1,0.00001,0.1))==90

end

@testset "Basic" begin



    
    function test_prepare_edges()
        ## Compared with pdelib; have more of these
        nodes=[
            0.0 1;
            1 0;
            0 1;
            1 1]'
        
        cells=[
            1 2 3;
            2 3 4
        ]'

        cellmat=[1,1]
        bfaces=[ 1 2;
                 1 3;
                 2 4]'

        bfacemat=[1,1]
        
        grid=simplexgrid(nodes,cells,cellmat,bfaces,bfacemat)

        @show grid[CellFaces]
        @show grid[FaceNodes]
        @show grid[FaceCells]
        
        @show grid[CellEdges]
        @show grid[EdgeNodes]
        @show grid[EdgeCells]
        
        grid[CellEdges]==Int32[1 2; 2 4; 3 5] &&       
        grid[FaceNodes]==Int32[1 2 3 3 4; 2 3 1 4 2] &&
            grid[EdgeCells]==Int32[1 1 1 2 2; 0 2 0 0 0]

    end
    @test test_prepare_edges()


    @test let
        g=simplexgrid(0:10)
        bfc=g[BFaceCells]
        bfc[1,1]==1 && bfc[1,2]==10
    end

    @test let
        g=simplexgrid(0:10)
        bfn=g[BFaceNormals]
        bfn==[-1.0 1.0]
    end

    
    @test let
        g=simplexgrid(0:2, 0:2)
        bfn=g[BFaceNormals]
        bfn==[0.0  -1.0;
              -1.0   0.0;
              0.0  -1.0;
              1.0   0.0;
              0.0   1.0;
              -1.0   0.0;
              1.0   0.0;
              0.0   1.0]'
    end

    @test let
        g=simplexgrid(0:2, 0:2)
        bfc=g[BFaceCells]
        tryfix(bfc)==[1  2  3  3  6  6  7  8;]
    end

    @test let
        X=[0 0; 2 0 ; 1 2;  1 0.5]'
        C=[ 1 2 4 ; 2 3 4 ; 3 1 4]'
        CR=ones(Int,3)
        F=[ 1 2 ; 2 3 ; 3 1]'
        FR=[1,2,3]
        g=simplexgrid(X,C,CR,F,FR)
        bfc=g[BFaceCells]
        tryfix(bfc)==[1  2  3;]
    end

    
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


@testset "1D" begin
    @test testgrid(interval_from_vector(),(21,20,2))
    @test testgrid(interval_localref(),   (27,26,2))
    @test testgrid(interval_multiregion(),(21,20,3))
    @test testgrid(interval_subgrid(),(51,50,2))
end

@testset "2D" begin
    @test testgrid(rectangle(),(441,800,80))
    @test testgrid(rectangle_localref(),(729, 1352, 104))
    @test testgrid(rectangle_multiregion(),(441,800,100))
    @test testgrid(rectangle_subgrid(),(360, 600, 120))
end

@testset "3D" begin
    @test testgrid(quadrilateral(),(330,1200,440))
    @test mask_bedges()
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


@testset "Grid Stuff" begin
    include("test_gridstuff.jl")
    run_grid_tests()
end

