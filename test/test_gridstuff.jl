function get_testgrid(::Type{<:Edge1D})
    X=collect(0:0.05:1)
    simplexgrid(X)
end

function get_testgrid(::Type{<:Triangle2D})
    X=collect(0:0.05:1)
    Y=collect(0:0.05:1)
    simplexgrid(X,Y)
end

function get_testgrid(::Type{<:Tetrahedron3D})
    X=collect(0:0.1:1)
    Y=collect(0:0.1:1)
    Z=collect(0:0.1:1)
    simplexgrid(X,Y,Z)
end

function get_testgrid(::Type{Parallelogram2D})
    return uniform_refine(grid_unitsquare(Parallelogram2D),2)
end

function get_testgrid(::Type{Parallelepiped3D})
    return uniform_refine(grid_unitcube(Parallelepiped3D),2)
end

function get_testgrid(::Type{Triangle2D},::Type{Parallelogram2D})
    return uniform_refine(grid_unitsquare_mixedgeometries(),1)
end

function check_enumeration_consistency(G::Type{<:AbstractElementGeometry}, CellItems, ItemNodes, item_enum_rule)

    @info "Checking enumeration consistency $CellItems <> $ItemNodes with $item_enum_rule for geometry = $G..."

    xgrid = get_testgrid(G)

    xCellNodes = xgrid[CellNodes]
    xCellItems = xgrid[CellItems]
    xItemNodes = xgrid[ItemNodes]
    item_rule = item_enum_rule(G)
    nitems_per_cell = max_num_targets_per_source(xCellItems)
    nnodes_per_item = max_num_targets_per_source(xItemNodes)
    ncells = num_sources(xCellNodes)

    @assert nitems_per_cell == size(item_rule,1)
    @assert nnodes_per_item == size(item_rule,2)

    item::Int = 0
    everything_is_consistent = true
    itemnodes_found = zeros(Bool, nnodes_per_item)
    for cell = 1 : ncells
        for i = 1 : nitems_per_cell
            item = xCellItems[i,cell]
            fill!(itemnodes_found,false)
            for n = 1 : nnodes_per_item
                for k = 1 : nnodes_per_item
                    if xItemNodes[n,item] == xCellNodes[item_rule[i,k],cell]
                        itemnodes_found[n] = true
                        break
                    end
                end
            end
            if prod(itemnodes_found) != true
                everything_is_consistent = false
            end
        end
    end

    return everything_is_consistent
end


function check_cellfinder(xgrid)
    EG = xgrid[UniqueCellGeometries]
    @info "Testing CellFinder for geometries=$EG..."
    xCoordinates = xgrid[Coordinates]
    xCellNodes = xgrid[CellNodes]
    edim = dim_element(EG[1])
    CF::CellFinder{Float64,Int32} = CellFinder(xgrid)
    cell_to_find = num_sources(xCellNodes) - 1

    # compute midpoint of cell_to_find
    x_source = zeros(Float64,edim)
    for j = 1 : edim, k = 1 : num_targets(xCellNodes,cell_to_find)
        x_source[j] += xCoordinates[j,xCellNodes[k,cell_to_find]] + 1e-6
    end
    x_source ./= num_targets(xCellNodes,cell_to_find)

    # find cell by local strategy
    xref = zeros(Float64,edim+1)
    cell = gFindLocal!(xref, CF, x_source; icellstart = 1)

    # check xref
    x = zeros(Float64,edim)
    L2G = L2GTransformer{Float64,xgrid[CellGeometries][cell],xgrid[CoordinateSystem]}(xgrid, ON_CELLS)
    update!(L2G,cell)
    eval!(x,L2G,xref)
    
    @info "... found x=$x in cell = $cell by local search (and had to find x=$x_source in cell=$cell_to_find)"

    @assert cell == cell_to_find

    # find cell again by brute force
    cell = gFindBruteForce!(xref, CF, x_source)

    # check xref
    x = zeros(Float64,edim)
    L2G = L2GTransformer{Float64,xgrid[CellGeometries][cell],xgrid[CoordinateSystem]}(xgrid, ON_CELLS)
    update!(L2G,cell)
    eval!(x,L2G,xref)
    
    @info "... found x=$x in cell = $cell by brute force (and had to find  x=$x_source in cell=$cell_to_find)"

    @assert cell == cell_to_find

    return (cell == cell_to_find) && (sqrt(sum((x - x_source).^2)) < 1e-14)
end


function run_grid_tests()
    # test FACE enumerations
    @test check_enumeration_consistency(Edge1D, CellFaces, FaceNodes, face_enum_rule)
    @test check_enumeration_consistency(Triangle2D, CellFaces, FaceNodes, face_enum_rule)
    @test check_enumeration_consistency(Parallelogram2D, CellFaces, FaceNodes, face_enum_rule)
    @test check_enumeration_consistency(Tetrahedron3D, CellFaces, FaceNodes, face_enum_rule)
    
    # test EDGE enumerations
    @test check_enumeration_consistency(Tetrahedron3D, CellEdges, EdgeNodes, edge_enum_rule)

    # todo: check FaceEdges, CellFaceSigns, CellFaceOrientations
    # todo: FaceCells, EdgeCells


    # @test check_cellfinder(get_testgrid(Edge1D))
    # @test check_cellfinder(get_testgrid(Triangle2D))
    # @test check_cellfinder(get_testgrid(Parallelogram2D))
    # @test check_cellfinder(get_testgrid(Triangle2D,Parallelogram2D))
    # @test check_cellfinder(get_testgrid(Tetrahedron3D))
    # @test check_cellfinder(get_testgrid(Parallelepiped3D))
    # @test check_cellfinder(get_testgrid(Triangle2D))

end
