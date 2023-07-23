ENV["MPLBACKEND"]="agg"

using Test, ExtendableGrids, GridVisualize, SHA, SimplexGridFactory, Triangulate
import Gmsh: gmsh # trigger extension
import CairoMakie

CairoMakie.activate!(type="svg",visible=false)


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

function testrw(grid,format;confidence=:medium)
    ftmp=tempname()*"."*format
    write(ftmp,grid)
    grid1=simplexgrid(ftmp)
    seemingly_equal(grid1,grid;confidence=confidence)
end


@testset "Read/Write sg" begin
    X=collect(0:0.05:1)
    @test testrw(simplexgrid(X),"sg")
    @test testrw(simplexgrid(X,X),"sg")
    @test testrw(simplexgrid(X,X,X),"sg")
end

@testset "Read/Write gmsh" begin
    X=collect(0:0.05:1)
    gmsh.initialize()
    @test testrw(simplexgrid(X,X),"msh";confidence=:low)
    @test testrw(simplexgrid(X,X,X),"msh";confidence=:low)
    gmsh.finalize()
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
    function rect1d()
        X=collect(0:0.1:1)
        g=simplexgrid(X)
        rect!(g, [0.3],[0.6], region=2,bregions=[3,4])
    end
    
    @test testgrid(interval_from_vector(),(21,20,2))
    @test testgrid(interval_localref(),   (27,26,2))
    @test testgrid(interval_multiregion(),(21,20,3))
    @test testgrid(interval_subgrid(),(51,50,2))
    @test testgrid(rect1d(),(11,10,4))
end

@testset "2D" begin
    function rect2d()
        X=collect(0:0.1:1)
        g=simplexgrid(X,X)
        rect!(g, [0.3, 0.3],[0.6,0.6], region=2,bregions=[3,4,5,6])
    end

    @test testgrid(rectangle(),(441,800,80))
    @test testgrid(rectangle_localref(),(729, 1352, 104))
    @test testgrid(rectangle_multiregion(),(441,800,100))
    @test testgrid(rectangle_subgrid(),(360, 600, 120))
    @test testgrid(rect2d(),(121,200,52))
    @test testgrid(rect2d_bregion_function(),(79,112,44))

    g,sg,sf=sorted_subgrid()
    @test testgrid(g,(187,306,66))
    @test testgrid(sg,(17,16,0))
    @test issorted(view(sg[Coordinates],1,:))

end

@testset "3D" begin
    function subgen(;h=0.2)
        X=collect(-1:h:2)
        Y=collect(0:h:1)
        Z=collect(0:h:2)
        g=simplexgrid(X,Y,Z)
        # Mark all boundaries as region 1
        bfacemask!(g,[-1,0,0], [2,1,2],1, allow_new=false)
        
        # Default cell region is 1
        # Mark cut-out regions as 2
        cellmask!(g, [-1,0,0], [0,1,1], 2)
        cellmask!(g, [1,0,0], [2,1,1], 2)
        
        # add new interface elements
        bfacemask!(g, [-1,0,1], [2,1,1],3,allow_new=true)
        bfacemask!(g, [0,0,0], [0,1,1],3,allow_new=true)
        bfacemask!(g, [1,0,0], [1,1,1],3,allow_new=true)
        
        # return subgrid of region 1
        subgrid(g,[1])
    end

    function rect3d()
        X=collect(0:0.1:1)
        g=simplexgrid(X,X,X)
        rect!(g, [0.3, 0.3,0.4],[0.6,0.6,0.7], region=2,bregion=8)
        g
    end


    
    @test testgrid(quadrilateral(),(330,1200,440))
    @test mask_bedges()

    X=collect(0:0.25:1)
    gxy=simplexgrid(X,X)
    gxyz=simplexgrid(gxy,X)
    g=simplexgrid(X,X,X)
    @test testgrid(gxyz,(125,384,192))
    @test g[Coordinates]≈gxyz[Coordinates]
    @test testgrid(subgen(),(756,3000,950))
    @test testgrid(rect3d(),(1331,6000,1308))
    @test testgrid(cross3d(),(189,480,344))
end

@testset "quadmeshes" begin
    function try_Quad()
        Tc = Float32
        Ti = Int32
        grid = ExtendableGrid{Tc, Ti}()
        
        grid[Coordinates] = convert(Matrix{Tc}, [0.0 0.0; 0.0 1.0; 1.0 0.0; 1.0 1.0]')
        grid[CellNodes]   = convert(Matrix{Ti}, [1 3 4 2]')
        
        grid[CellRegions] = convert(Vector{Ti}, [1])
    
        
        grid[CellGeometries] = VectorOfConstants{ElementGeometries,Ti}(Quadrilateral2D,length(grid[CellRegions]))
        
        grid[BFaceNodes]  = convert(Matrix{Ti}, [1 3; 3 4; 4 2; 2 1]')
        grid[BFaceRegions] = convert(Vector{Ti}, [1, 1, 2, 2])
        
        grid
    end
    @test isa(try_Quad(), ExtendableGrid)
end

@testset "plotting examples" begin
    include("../docs/makeplots.jl")
    picdir=mktempdir()
    
    @test makeplot("interval_from_vector",picdir; Plotter=CairoMakie)
    @test makeplot("interval_localref",picdir; Plotter=CairoMakie)
    @test makeplot("interval_multiregion",picdir; Plotter=CairoMakie)
    @test makeplot("interval_subgrid",picdir; Plotter=CairoMakie)
    @test makeplot("rectangle",picdir; Plotter=CairoMakie)
    @test makeplot("rectangle_localref",picdir; Plotter=CairoMakie)
    @test makeplot("rectangle_multiregion",picdir; Plotter=CairoMakie)
    @test makeplot("rectangle_subgrid",picdir; Plotter=CairoMakie)
    @test makeplot("quadrilateral",picdir; Plotter=CairoMakie)        
    @test makeplotx("sorted_subgrid",picdir; Plotter=CairoMakie)        
end



function tglue(;dim=2,breg=0)
    X1=linspace(0,1,5)
    Y1=linspace(0,1,5)
    Z1=linspace(0,1,5)

    X2=linspace(1,2,5)
    Y2=linspace(0,0.5,3)
    Z2=linspace(0,0.5,3)

    if dim==1
        g1=simplexgrid(X1)
        g2=simplexgrid(X2)
    end
    
    if dim==2
        g1=simplexgrid(X1,Y1)
        g2=simplexgrid(X2,Y2)
    end
    
    if dim==3
        g1=simplexgrid(X1,Y1,Z1)
        g2=simplexgrid(X2,Y2,Z2)
    end
    
    glue(g1,g2,interface=breg)
end

@testset "Glue" begin
    @test testgrid(tglue(;dim=1,breg=0),(9,8,2))
    @test testgrid(tglue(;dim=1,breg=1),(9,8,3))
    @test testgrid(tglue(;dim=2,breg=0),(37, 48, 24))
    @test testgrid(tglue(;dim=2,breg=1),(37, 48, 26))
    @test testgrid(tglue(;dim=3,breg=0),(161, 480, 256))
    @test testgrid(tglue(;dim=3,breg=1),(161, 480, 264))
end


    
function voronoitest()
    g=simplexgrid(0:0.1:1)
    @test g[VoronoiFaceCenters]≈[0.05  0.15  0.25  0.35  0.45  0.55  0.65  0.75  0.85  0.95]
    g=simplexgrid(0:0.5:1,0:0.5:1)
    @test g[VoronoiFaceCenters]≈[0.25  0.5   0.25  0.25  0.0   0.75  1.0   0.75  0.75  0.5   0.25  0.25  0.0   1.0   0.75  0.75;
                                  0.0   0.25  0.25  0.5   0.25  0.0   0.25  0.25  0.5   0.75  0.75  1.0   0.75  0.75  0.75  1.0]
end

@testset "Voronoi" begin
    voronoitest()
end
                           
@testset "Grid Stuff" begin
    include("test_gridstuff.jl")
    run_grid_tests()
end         

@testset "Interpolations" begin
    X_from = collect(0:0.1:1)
    X_to = collect(0:0.02:1)

    grid_from=simplexgrid(X_from)
    grid_to=simplexgrid(X_to)
    u_from=map((x)->x,grid_from)
    u_to_exact=map((x)->x,grid_to)
    u_to=interpolate(grid_to,u_from,grid_from)
    @test u_to≈u_to_exact
    
    u_from=Matrix(hcat(u_from,u_from)')
    u_to_exact=Matrix(hcat(u_to_exact,u_to_exact)')
    u_to=interpolate(grid_to,u_from,grid_from)
    @test u_to≈u_to_exact
    
    grid_from=simplexgrid(X_from,X_from)
    grid_to=simplexgrid(X_to,X_to)
    u_from=map((x,y)->x+y,grid_from)
    u_to_exact=map((x,y)->x+y,grid_to)
    u_to=interpolate(grid_to,u_from,grid_from)
    @test u_to≈u_to_exact
    
    u_from=Matrix(hcat(u_from,u_from)')
    u_to_exact=Matrix(hcat(u_to_exact,u_to_exact)')
    u_to=interpolate(grid_to,u_from,grid_from)
    @test u_to≈u_to_exact
    
    grid_from=simplexgrid(X_from,X_from,X_from)
    grid_to=simplexgrid(X_to,X_to,X_to)
    u_from=map((x,y,z)->x+y+z,grid_from)
    u_to_exact=map((x,y,z)->x+y+z,grid_to)
    u_to=interpolate(grid_to,u_from,grid_from)
    @test u_to≈u_to_exact
    
    u_from=Matrix(hcat(u_from,u_from)')
    u_to_exact=Matrix(hcat(u_to_exact,u_to_exact)')
    u_to=interpolate(grid_to,u_from,grid_from)
    @test u_to≈u_to_exact
    
end

@testset "writeVTK" begin
    X = collect(-1:1.0:1)
    Y = collect(-1:1.0:1)
    Z = collect(-2:1.0:2)
    g = simplexgrid(X,Y,Z)

    nx = num_nodes(g)
    nc = num_cells(g)

    # ensure calculation of these data is free of roundoff errors
    point_data    = map((x,y,z) -> (x+y+z), g)
    field_data = [1.0, 2, 3, 4, 5, 6]

    writeVTK("testfile_writevtk.vtu", g; 
        cellregions = g[CellRegions],
        point_data = point_data, 
        field_data = field_data)

    sha_code = open("testfile_writevtk.vtu") do f
        sha256(f)
    end |> bytes2hex

    @test sha_code == "93a31139ccb3ae3017351d7cef0c2639c5def97c9744699543fe8bc58e1ebcea"
end

