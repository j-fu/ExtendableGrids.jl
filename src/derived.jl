# Here goes all the stuff that connects this code to ExtendableGrids like
# - definition of additional adjacency types and their instantiation
# - connections between CellGeometries and FaceGeometries of their faces
# - formulas to compute Volumes and Normal Vectors etc.
# (some of this might become native in the ExtendableGrids module itself at some point)


# additional ElementGeometryTypes with parent information
#abstract type Vertex0DWithParent{Parent <: AbstractElementGeometry} <: Vertex0D end
#abstract type Vertex0DWithParents{Parent1 <: AbstractElementGeometry, Parent2 <: AbstractElementGeometry} <: Vertex0D end
#export Vertex0DWithParent, Vertex0DWithParents

#function AddParent(FEG::Type{<:Vertex0D}, CEG::Type{<:AbstractElementGeometry})
#    return Vertex0DWithParent{CEG}
#end

#abstract type Edge1DWithParent{Parent <: AbstractElementGeometry} <: Edge1D end
#abstract type Edge1DWithParents{Parent1 <: AbstractElementGeometry, Parent2 <: AbstractElementGeometry} <: Edge1D end
#export Edge1DWithParent, Edge1DWithParents

#function AddParent(FEG::Type{<:Edge1D}, CEG::Type{<:AbstractElementGeometry})
#    return Edge1DWithParent{CEG}
#end

const GridEGTypes = Vector{DataType}


# additional ExtendableGrids adjacency types 
abstract type CellEdges <: AbstractGridAdjacency end
abstract type CellFaces <: AbstractGridAdjacency end
abstract type CellFaceSigns <: AbstractGridAdjacency end
abstract type CellFaceOrientations <: AbstractGridAdjacency end
abstract type CellEdgeSigns <: AbstractGridAdjacency end
abstract type CellVolumes <: AbstractGridFloatArray1D end
abstract type UniqueCellGeometries <: AbstractElementGeometries end

abstract type FaceNodes <: AbstractGridAdjacency end
abstract type FaceVolumes <: AbstractGridFloatArray1D end
abstract type FaceCells <: AbstractGridAdjacency end
abstract type FaceEdges <: AbstractGridAdjacency end
abstract type FaceNormals <: AbstractGridFloatArray2D end
abstract type FaceGeometries <: AbstractElementGeometries end
abstract type FaceRegions <: AbstractElementRegions end
abstract type UniqueFaceGeometries <: AbstractElementGeometries end

abstract type BFaces <: AbstractGridIntegerArray1D end
abstract type BFaceCellPos <: AbstractGridIntegerArray1D end # position of bface in adjacent cell
abstract type BFaceVolumes <: AbstractGridFloatArray1D end
abstract type UniqueBFaceGeometries <: AbstractElementGeometries end

abstract type EdgeNodes <: AbstractGridAdjacency end
abstract type EdgeVolumes <: AbstractGridFloatArray1D end
abstract type EdgeCells <: AbstractGridAdjacency end
abstract type EdgeTangents <: AbstractGridFloatArray2D end
abstract type EdgeRegions <: AbstractElementRegions end
abstract type EdgeGeometries <: AbstractElementGeometries end
abstract type UniqueEdgeGeometries <: AbstractElementGeometries end

abstract type BEdgeNodes <: AbstractGridAdjacency end
abstract type BEdges <: AbstractGridIntegerArray1D end
#abstract type BEdgeCellPos <: AbstractGridIntegerArray1D end # position of bface in adjacent cell
abstract type BEdgeVolumes <: AbstractGridFloatArray1D end
abstract type BEdgeRegions <: AbstractElementRegions end
abstract type BEdgeGeometries <: AbstractElementGeometries end
abstract type UniqueBEdgeGeometries <: AbstractElementGeometries end

abstract type NodePatchGroups <: AbstractGridIntegerArray1D end


## grid item types to dispatch certain things to the correct GridComponents
abstract type ITEMTYPE_NODE end
abstract type ITEMTYPE_CELL end
abstract type ITEMTYPE_FACE end
abstract type ITEMTYPE_BFACE end
abstract type ITEMTYPE_EDGE end
abstract type ITEMTYPE_BEDGE end

abstract type PROPERTY_NODES end
abstract type PROPERTY_VOLUME end
abstract type PROPERTY_REGION end
abstract type PROPERTY_GEOMETRY end
abstract type PROPERTY_UNIQUEGEOMETRY end

GridComponent4TypeProperty(::Type{ITEMTYPE_CELL},::Type{PROPERTY_NODES}) = CellNodes
GridComponent4TypeProperty(::Type{ITEMTYPE_CELL},::Type{PROPERTY_VOLUME}) = CellVolumes
GridComponent4TypeProperty(::Type{ITEMTYPE_CELL},::Type{PROPERTY_REGION}) = CellRegions
GridComponent4TypeProperty(::Type{ITEMTYPE_CELL},::Type{PROPERTY_GEOMETRY}) = CellGeometries
GridComponent4TypeProperty(::Type{ITEMTYPE_CELL},::Type{PROPERTY_UNIQUEGEOMETRY}) = UniqueCellGeometries

GridComponent4TypeProperty(::Type{ITEMTYPE_FACE},::Type{PROPERTY_NODES}) = FaceNodes
GridComponent4TypeProperty(::Type{ITEMTYPE_FACE},::Type{PROPERTY_VOLUME}) = FaceVolumes
GridComponent4TypeProperty(::Type{ITEMTYPE_FACE},::Type{PROPERTY_REGION}) = FaceRegions
GridComponent4TypeProperty(::Type{ITEMTYPE_FACE},::Type{PROPERTY_GEOMETRY}) = FaceGeometries
GridComponent4TypeProperty(::Type{ITEMTYPE_FACE},::Type{PROPERTY_UNIQUEGEOMETRY}) = UniqueFaceGeometries

GridComponent4TypeProperty(::Type{ITEMTYPE_BFACE},::Type{PROPERTY_NODES}) = BFaceNodes
GridComponent4TypeProperty(::Type{ITEMTYPE_BFACE},::Type{PROPERTY_VOLUME}) = BFaceVolumes
GridComponent4TypeProperty(::Type{ITEMTYPE_BFACE},::Type{PROPERTY_REGION}) = BFaceRegions
GridComponent4TypeProperty(::Type{ITEMTYPE_BFACE},::Type{PROPERTY_GEOMETRY}) = BFaceGeometries
GridComponent4TypeProperty(::Type{ITEMTYPE_BFACE},::Type{PROPERTY_UNIQUEGEOMETRY}) = UniqueBFaceGeometries

GridComponent4TypeProperty(::Type{ITEMTYPE_EDGE},::Type{PROPERTY_NODES}) = EdgeNodes
GridComponent4TypeProperty(::Type{ITEMTYPE_EDGE},::Type{PROPERTY_VOLUME}) = EdgeVolumes
GridComponent4TypeProperty(::Type{ITEMTYPE_EDGE},::Type{PROPERTY_REGION}) = EdgeRegions
GridComponent4TypeProperty(::Type{ITEMTYPE_EDGE},::Type{PROPERTY_GEOMETRY}) = EdgeGeometries
GridComponent4TypeProperty(::Type{ITEMTYPE_EDGE},::Type{PROPERTY_UNIQUEGEOMETRY}) = UniqueEdgeGeometries

GridComponent4TypeProperty(::Type{ITEMTYPE_BEDGE},::Type{PROPERTY_NODES}) = BEdgeNodes
GridComponent4TypeProperty(::Type{ITEMTYPE_BEDGE},::Type{PROPERTY_VOLUME}) = BEdgeVolumes
GridComponent4TypeProperty(::Type{ITEMTYPE_BEDGE},::Type{PROPERTY_REGION}) = BEdgeRegions
GridComponent4TypeProperty(::Type{ITEMTYPE_BEDGE},::Type{PROPERTY_GEOMETRY}) = BEdgeGeometries
GridComponent4TypeProperty(::Type{ITEMTYPE_BEDGE},::Type{PROPERTY_UNIQUEGEOMETRY}) = UniqueBEdgeGeometries


function get_facegrid(source_grid::ExtendableGrid{Tc,Ti}) where {Tc,Ti}
    facegrid=ExtendableGrid{Tc,Ti}()
    facegrid[Coordinates]=source_grid[Coordinates]
    facegrid[CellNodes]=source_grid[FaceNodes]
    facegrid[CoordinateSystem]=source_grid[CoordinateSystem]
    facegrid[CellGeometries]=source_grid[FaceGeometries]
    facegrid[UniqueCellGeometries]=source_grid[UniqueFaceGeometries]
    facegrid[CellRegions] = source_grid[FaceRegions]
    # todo: facegrid[CellFaces] = source_grid[FaceEdges]
    return facegrid
end

function get_bfacegrid(source_grid::ExtendableGrid{Tc,Ti})  where {Tc,Ti}
    bfacegrid=ExtendableGrid{Tc,Ti}()
    bfacegrid[Coordinates]=source_grid[Coordinates]
    bfacegrid[CellNodes]=source_grid[BFaceNodes]
    bfacegrid[CoordinateSystem]=source_grid[CoordinateSystem]
    bfacegrid[CellGeometries]=source_grid[BFaceGeometries]
    bfacegrid[UniqueCellGeometries]=source_grid[UniqueBFaceGeometries]
    bfacegrid[CellRegions] = source_grid[BFaceRegions]
    # todo: bfacegrid[CellFaces] = source_grid[BFaceEdges] (or sub-view of FaceEdges ?)
    return bfacegrid
end

function get_edgegrid(source_grid::ExtendableGrid{Tc,Ti}) where {Tc,Ti}
    edgegrid=ExtendableGrid{Tc,Ti}()
    edgegrid[Coordinates]=source_grid[Coordinates]
    edgegrid[CellNodes]=source_grid[EdgeNodes]
    edgegrid[CoordinateSystem]=source_grid[CoordinateSystem]
    edgegrid[CellGeometries]=source_grid[EdgeGeometries]
    edgegrid[UniqueCellGeometries]=source_grid[UniqueedgeGeometries]
    edgegrid[CellRegions] = source_grid[EdgeRegions]
    return edgegrid
end

# show function for ExtendableGrids and defined Components in its Dict
function showmore(io::IO, xgrid::ExtendableGrid{Tc,Ti}) where {Tc,Ti}

    dim = size(xgrid[Coordinates],1)
    nnodes = num_sources(xgrid[Coordinates])
    ncells = num_sources(xgrid[CellNodes])
    
	println("ExtendableGrid information");
    println("==========================");
	println("dim: $(dim)")
	println("nnodes: $(nnodes)")
    println("ncells: $(ncells)")
    if haskey(xgrid.components,FaceNodes)
        nfaces = num_sources(xgrid[FaceNodes])
        println("nfaces: $(nfaces)")
    else
        println("nfaces: (FaceNodes not instantiated)")
    end
    if haskey(xgrid.components,EdgeNodes)
        nfaces = num_sources(xgrid[FaceNodes])
        println("nedges: $(nedges)")
    else
        println("nedges: (EdgeNodes not instantiated)")
    end
    println("")
    println("Components");
    println("==========");
    for tuple in xgrid.components
        println("> $(tuple[1])")
    end
end


# FaceNodes = nodes for each face (implicitly defines the enumerations of faces)
function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{FaceNodes}) where {Tc,Ti}

    xCellNodes = xgrid[CellNodes]
    ncells = num_sources(xCellNodes)
    nnodes = num_sources(xgrid[Coordinates])
    xCellGeometries = xgrid[CellGeometries]

    # transpose CellNodes to get NodeCells
    xNodeCells = atranspose(xCellNodes)
    max_ncell4node::Ti = max_num_targets_per_source(xNodeCells)

    # instantiate new empty adjacency fields
    xFaceCells = zeros(Ti,0) # cells are appended and at the end rewritten into 2,nfaces array


    # find unique face enumeration rules
    EG = unique(xCellGeometries)
    dim::Ti = dim_element(EG[1])
    face_rules = Array{Array{Ti,2},1}(undef,length(EG))
    maxfacenodes::Ti = 0
    FEG = []
    for j = 1 : length(EG)
        face_rules[j] = face_enum_rule(EG[j])
        maxfacenodes = max(size(face_rules[j],2),maxfacenodes)
        for k = 1 : nfaces_for_geometry(EG[j])
            append!(FEG,[facetype_of_cellface(EG[j],j)])
        end
    end
    unique!(FEG)

    # check if only one type of cell geometry is present in grid
    if length(EG) == 1
        singleEG = true 
    else
        singleEG = false
    end
    # check if only one type of face geometry is present in grid
    if length(FEG) == 1
        singleFEG = true
    else
        singleFEG = false
    end
    
    if singleFEG
        xFaceNodes = zeros(Ti,0)
    else
        xFaceNodes = VariableTargetAdjacency(Ti)
        xFaceGeometries::Array{DataType,1} = []
    end
    if singleEG == true && singleFEG == true
        # only one geometry type allows for much faster code
        xCellFaces = zeros(Ti,nfaces_for_geometry(EG[1]),ncells)
        xCellFaceSigns = zeros(Ti,nfaces_for_geometry(EG[1]),ncells)
    else
        xCellFaces = VariableTargetAdjacency(Ti)
        xCellFaceSigns = VariableTargetAdjacency(Ti)
        # pre-allocate xCellFaces
        cellEG = xCellGeometries[1]
        for cell = 1 : ncells
            cellEG = xCellGeometries[cell]
            append!(xCellFaces,zeros(Ti,nfaces_for_geometry(cellEG)))
            append!(xCellFaceSigns,zeros(Ti,nfaces_for_geometry(cellEG)))
        end   
    end

    # temporary variables
    # pre-initialised ready to work for singleEG and singleFEG
    node::Ti = 0
    face::Ti = 0
    cell::Ti = 0
    cell2::Ti = 0
    cellEG = EG[1]
    cell2EG = EG[1]
    faceEG = FEG[1]
    faceEG2 = FEG[1]
    face_rule::Array{Ti,2} = face_rules[1]
    face_rule2::Array{Ti,2} = face_rules[1]
    iEG::Ti = 1
    nneighbours::Ti = 0
    faces_per_cell::Ti = nfaces_for_geometry(cellEG)
    faces_per_cell2::Ti = nfaces_for_geometry(cellEG)
    nodes_per_cellface::Ti = nnodes_for_geometry(faceEG)
    current_item::Array{Ti,1} = zeros(Ti,maxfacenodes) # should be large enough to store largest nnodes per cellface
    flag4item::Array{Bool,1} = zeros(Bool,nnodes)
    no_neighbours_found::Bool = true
    same_face::Bool = false
    
    # loop over cells
    for cell = 1 : ncells
        # find EG index for geometry
        if !singleEG
            cellEG = xCellGeometries[cell]
            faces_per_cell = nfaces_for_geometry(cellEG)
            for j=1:length(EG)
                if cellEG == EG[j]
                    iEG = j
                    break;
                end
            end
            face_rule = face_rules[iEG]
        end

        # loop over cell faces
        for k = 1 : faces_per_cell

            # check if face is already known to cell
            if xCellFaces[k,cell] > 0
                continue;
            end    

            # get face geometry
            if !singleFEG
                faceEG = facetype_of_cellface(cellEG, k)
                nodes_per_cellface = nnodes_for_geometry(faceEG)
            end

            # flag face nodes and commons4cells
            for j = nodes_per_cellface:-1:1
                node = xCellNodes[face_rule[k,j],cell]
                current_item[j] = node
                flag4item[node] = true; 
            end

            # get neighbours for first node
            nneighbours = num_targets(xNodeCells,node)

            # loop over neighbours
            no_neighbours_found = true
            for n = 1 : nneighbours
                cell2 = xNodeCells[n,node]

                # skip if cell2 is the same as cell
                if (cell == cell2) 
                    continue; 
                end

                # find face enumeration rule
                if !singleEG
                    cell2EG = xCellGeometries[cell2]
                    faces_per_cell2 = nfaces_for_geometry(cell2EG)
                    for j=1:length(EG)
                        if cell2EG == EG[j]
                            iEG = j
                            break;
                        end
                    end
                    face_rule2 = face_rules[iEG]
                end

                # loop over faces face2 of adjacent cell2
                for f2 = 1 : faces_per_cell2

                    # check if face f2 is already known to cell2
                    if xCellFaces[f2,cell2] != 0
                        continue;
                    end    

                    # check if face f2 has same geometry
                    if !singleFEG
                        faceEG2 = facetype_of_cellface(cell2EG, f2)
                        if faceEG != faceEG2
                            continue;
                        end
                    end

                    # otherwise compare nodes of face and face2
                    same_face = true
                    for j = 1 : nodes_per_cellface
                        if flag4item[xCellNodes[face_rule2[f2,j],cell2]] == false
                            same_face = false
                            break;
                        end    
                    end
                    
                    # if all nodes are the same, register face
                    if (same_face)
                        no_neighbours_found = false
                        face += 1
                        push!(xFaceCells,cell)
                        push!(xFaceCells,cell2)
                        if singleEG == false
                            xCellFaces.colentries[xCellFaces.colstart[cell]+k-1] = face
                            xCellFaces.colentries[xCellFaces.colstart[cell2]+f2-1] = face
                            xCellFaceSigns.colentries[xCellFaceSigns.colstart[cell]+k-1] = 1
                            xCellFaceSigns.colentries[xCellFaceSigns.colstart[cell2]+f2-1] = -1
                            if singleFEG == false
                                push!(xFaceGeometries,faceEG)
                            end
                        else
                            xCellFaces[k,cell] = face
                            xCellFaces[f2,cell2] = face
                            xCellFaceSigns[k,cell] = 1
                            xCellFaceSigns[f2,cell2] = -1
                        end
                        append!(xFaceNodes,view(current_item,1:nodes_per_cellface))
                        break;
                    end
                end
            end

            # if no common neighbour cell is found, register face (boundary faces)
            if no_neighbours_found == true
                face += 1
                push!(xFaceCells,cell)
                push!(xFaceCells,0)
                if singleEG == false
                    xCellFaces.colentries[xCellFaces.colstart[cell]+k-1] = face
                    xCellFaceSigns.colentries[xCellFaceSigns.colstart[cell]+k-1] = 1
                    if singleFEG == false
                        push!(xFaceGeometries,faceEG)
                    end
                else
                    xCellFaces[k,cell] = face
                    xCellFaceSigns[k,cell] = 1
                end
                append!(xFaceNodes,view(current_item,1:nodes_per_cellface))
            end

            # reset flag4item
            for j = 1:nodes_per_cellface
                flag4item[current_item[j]] = false 
            end
        end    
    end

    if singleFEG
        xFaceNodes = reshape(xFaceNodes,nodes_per_cellface,Ti(length(xFaceNodes)/nodes_per_cellface))
        xgrid[FaceGeometries] = VectorOfConstants(facetype_of_cellface(EG[1], 1), face)
    else
        xgrid[FaceGeometries] = xFaceGeometries
    end
    xgrid[CellFaces] = xCellFaces
    xgrid[CellFaceSigns] = xCellFaceSigns
    xgrid[FaceCells] = reshape(xFaceCells,2,face)
    xFaceNodes
end


function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{NodePatchGroups}) where {Tc,Ti}
    xCellNodes::GridAdjacencyTypes = xgrid[CellNodes]
    xNodeCells = atranspose(xCellNodes)
    nnodes = size(xgrid[Coordinates],2)
    ncells = num_sources(xCellNodes)
    ncells4node::Ti = 0
    group4node = zeros(Ti,nnodes)
    cgroup::Ti = 0
    cell_in_group = zeros(Bool,ncells)
    take_into = false
    cgroup = 0
    while minimum(group4node) == 0
        cell_in_group .= false
        cgroup += 1
        for node = 1 : nnodes
            if group4node[node] == 0
                ncells4node = num_targets(xNodeCells,node)
                take_into = true
                for c = 1 : ncells4node
                    if cell_in_group[xNodeCells[c,node]] == true
                        take_into = false
                        break
                    end
                end
                if take_into
                    group4node[node] = cgroup
                    for c = 1 : ncells4node
                        cell_in_group[xNodeCells[c,node]] = true
                    end
                end
            end
        end
    end
    group4node
end



# FaceNodes = nodes for each face (implicitly defines the enumerations of faces)
function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{EdgeNodes}) where {Tc,Ti}

    xCellNodes = xgrid[CellNodes]
    ncells::Ti = num_sources(xCellNodes)
    nnodes::Ti = num_sources(xgrid[Coordinates])
    xCellGeometries = xgrid[CellGeometries]
    dim::Ti = dim_element(xCellGeometries[1])
    if dim == 1
        # do nothing in 1D
        return
    end

    if dim == 2

        # xgrid[CellEdges] = xgrid[CellFaces]
        # xgrid[EdgeCells] = xgrid[FaceCells]
        # xgrid[EdgeNodes] = xgrid[FaceNodes]
        return  prepare_edges!(xgrid) # keep compatible to VoronoiFVM
    end

    
    #transpose CellNodes to get NodeCells
    xNodeCells = atranspose(xCellNodes)
    max_ncell4node::Ti = max_num_targets_per_source(xNodeCells)

    xEdgeCells = VariableTargetAdjacency(Ti)


    # find unique edge enumeration rules
    EG::Array{DataType,1} = unique(xCellGeometries)
    edge_rules::Array{Array{Ti,2},1} = Array{Array{Ti,2},1}(undef,length(EG))
    EEG = [edgetype_of_celledge(EG[1], 1)]
    maxedgenodes::Ti = 0
    for j = 1 : length(EG)
        edge_rules[j] = edge_enum_rule(EG[j])
        maxedgenodes = max(size(edge_rules[j],2),maxedgenodes)
    end
    if length(EG) == 1
        singleEG = true
        singleEEG = true # even in 3D every edge is Edge1D
        xCellEdges = zeros(Ti,nedges_for_geometry(EG[1]),ncells)
        xCellEdgeSigns = zeros(Ti,nedges_for_geometry(EG[1]),ncells)
    else
        singleEG = false
        singleEEG = true # even in 3D every edge is Edge1D
        xCellEdges = VariableTargetAdjacency(Ti)
        # pre-allocate xCellEdges
        cellEG = xCellGeometries[1]
        for cell = 1 : ncells
            cellEG = xCellGeometries[cell]
            append!(xCellEdges,zeros(Ti,nedges_for_geometry(cellEG)))
            append!(xCellEdgeSigns,zeros(Ti,nedges_for_geometry(cellEG)))
        end   
    end
    if !singleEEG
        xEdgeNodes = VariableTargetAdjacency(Ti)
        xEdgeGeometries::Array{DataType,1} = []
    else
        xEdgeNodes = zeros(Ti,0)
    end

    edge_rule::Array{Ti,2} = edge_rules[1]
    edge_rule2::Array{Ti,2} = edge_rules[1]
    current_item::Array{Ti,1} = zeros(Ti,nnodes_for_geometry(EEG[1])) # should be large enough to store largest nnodes_per_celledge
    flag4item::Array{Bool,1} = zeros(Bool,nnodes)
    cellEG = EG[1]
    cell2EG = EG[1]
    node::Ti = 0
    node_cells::Array{Ti,1} = zeros(Ti,max_ncell4node) # should be large enough to store largest nnodes_per_celledge
    cell2::Ti = 0
    nneighbours::Ti = 0
    edges_per_cell::Ti = nedges_for_geometry(EG[1])
    edges_per_cell2::Ti = nedges_for_geometry(EG[1])
    nodes_per_celledge::Ti = nnodes_for_geometry(EEG[1])
    nodes_per_celledge2::Ti = nnodes_for_geometry(EEG[1])
    common_nodes::Ti = 0
    cells_with_common_edge::Array{Ti,1} = zeros(Ti,max_ncell4node)
    pos_in_cells_with_common_edge::Array{Ti,1} = zeros(Ti,max_ncell4node)
    sign_in_cells_with_common_edge::Array{Ti,1} = zeros(Ti,max_ncell4node)
    ncells_with_common_edge::Ti = 0
    edge::Ti = 0
    iEG::Ti = 0

    # loop over cells
    for cell = 1 : ncells

        # find EG index for geometry
        if singleEG == false
            cellEG = xCellGeometries[cell]
            edges_per_cell = nedges_for_geometry(cellEG)
            for j=1:length(EG)
                if cellEG == EG[j]
                    iEG = j
                    break;
                end
            end
            edge_rule = edge_rules[iEG]
        end

        # loop over cell edges
        for k = 1 : edges_per_cell

            # check if edge is already known to cell
            if xCellEdges[k,cell] > 0
                continue;
            end    
            #nodes_per_celledge = nnodes_for_geometry(edgetype_of_celledge(cellEG, k))
            ncells_with_common_edge = 1
            cells_with_common_edge[1] = cell
            pos_in_cells_with_common_edge[1] = k
            sign_in_cells_with_common_edge[1] = 1

            # flag edge nodes and commons4cells
            for j = 1 : nodes_per_celledge
                node = xCellNodes[edge_rule[k,j],cell]
                current_item[j] = node
                flag4item[node] = true; 
            end

            # get first node and its neighbours
            node = xCellNodes[edge_rule[k,1],cell]
            nneighbours = num_targets(xNodeCells,node)
            node_cells[1:nneighbours] = xNodeCells[:,node]

            # loop over neighbours
            for n = 1 : nneighbours
                cell2 = node_cells[n]

                # skip if cell2 is the same as cell
                if (cell == cell2) 
                    continue; 
                end

                # loop over edges of cell2
                if singleEG == false
                    cell2EG = xCellGeometries[cell2]
                    edges_per_cell2 = nedges_for_geometry(cell2EG)

                    # find edge enumeration rule
                    for j=1:length(EG)
                        if cell2EG == EG[j]
                            iEG = j
                            break;
                        end
                    end
                    edge_rule2 = edge_rules[iEG]
                end

                for f2 = 1 : edges_per_cell2
                    # compare nodes of edge and edge2
                    #nodes_per_celledge2 = nnodes_for_geometry(edgetype_of_celledge(cell2EG, f2))
                    common_nodes = 0
                    if nodes_per_celledge == nodes_per_celledge2
                        for j = 1 : nodes_per_celledge2
                            if flag4item[xCellNodes[edge_rule2[f2,j],cell2]]
                                common_nodes += 1
                            else
                                continue;    
                            end    
                        end
                    end

                    # if all nodes are the same, register edge
                    if (common_nodes == nodes_per_celledge2)
                        ncells_with_common_edge += 1
                        cells_with_common_edge[ncells_with_common_edge] = cell2
                        pos_in_cells_with_common_edge[ncells_with_common_edge] = f2
                        if xCellNodes[edge_rule2[f2,1],cell2] == current_item[1]
                            sign_in_cells_with_common_edge[ncells_with_common_edge] = 1
                        else
                            sign_in_cells_with_common_edge[ncells_with_common_edge] = -1
                        end
                    end
                end
            end

            # register edge
            edge += 1
            if singleEG
                for c = 1 : ncells_with_common_edge
                    xCellEdges[pos_in_cells_with_common_edge[c],cells_with_common_edge[c]] = edge
                    xCellEdgeSigns[pos_in_cells_with_common_edge[c],cells_with_common_edge[c]] = sign_in_cells_with_common_edge[c]
                end
            else
                for c = 1 : ncells_with_common_edge
                    xCellEdges.colentries[xCellEdges.colstart[cells_with_common_edge[c]]+pos_in_cells_with_common_edge[c]-1] = edge
                    xCellEdgeSigns.colentries[xCellEdgeSigns.colstart[cells_with_common_edge[c]]+pos_in_cells_with_common_edge[c],cells_with_common_edge[c]] = sign_in_cells_with_common_edge[c]
                end
            end
            append!(xEdgeCells,view(cells_with_common_edge,1:ncells_with_common_edge))
            if singleEEG == false
                push!(xEdgeGeometries,edgetype_of_celledge(cellEG,k))
                append!(xEdgeNodes,view(current_item,1:nodes_per_celledge))
            else
                append!(xEdgeNodes,current_item)
            end

            #reset flag4item
            for j = 1 : nodes_per_celledge
                flag4item[current_item[j]] = false; 
            end
        end    
    end
    xgrid[CellEdges] = xCellEdges
    xgrid[EdgeCells] = xEdgeCells
    xgrid[CellEdgeSigns] = xCellEdgeSigns
    if singleEEG
        xgrid[EdgeGeometries] = VectorOfConstants(edgetype_of_celledge(EG[1],1),edge)
        reshape(xEdgeNodes,nnodes_for_geometry(edgetype_of_celledge(EG[1],1)),edge)
    else
        xgrid[EdgeGeometries] = xEdgeGeometries
        xEdgeNodes
    end
end


function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{CellFaceOrientations}) where {Tc,Ti}
    xCellFaceSigns = xgrid[CellFaceSigns]
    xCellGeometries = xgrid[CellGeometries]
    xCellFaces = xgrid[CellFaces]
    xFaceNodes = xgrid[FaceNodes]
    xCellNodes = xgrid[CellNodes]
    ncells = num_sources(xCellNodes)

    EG = unique(xCellGeometries)
    face_rules = Array{Array{Ti,2},1}(undef,length(EG))
    maxfacenodes = 0
    for j = 1 : length(EG)
        face_rules[j] = face_enum_rule(EG[j])
        maxfacenodes = max(size(face_rules[j],2),maxfacenodes)
    end
    singleEG = false
    if typeof(xCellFaceSigns) <: VariableTargetAdjacency
        xCellFaceOrientations = deepcopy(xCellFaceSigns)
    else
        singleEG = true
        xCellFaceOrientations = zeros(Ti, size(xCellFaceSigns,1), size(xCellFaceSigns,2))
    end

    cellEG = EG[1]
    face_rule = face_rules[1]
    ncellfaces::Ti = 0
    nfacenodes::Ti = 0
    face::Ti = 0
    facenodes = zeros(Ti,maxfacenodes)
    found_configuration::Bool = false
    n::Ti = 0
    iEG::Ti = 0
    for cell = 1 : ncells
        cellEG = xCellGeometries[cell]

        # find EG index for geometry
        for j=1:length(EG)
            if cellEG == EG[j]
                iEG = j
                break;
            end
        end
        face_rule = face_rules[iEG] # determines local enumeration of faces

        # determine orientation
        ncellfaces = num_targets(xCellFaces,cell)
        for j = 1 : ncellfaces
            face = xCellFaces[j,cell]
            nfacenodes = num_targets(xFaceNodes,face)
            if xCellFaceSigns[j,cell] == 1
                if singleEG == false
                    xCellFaceOrientations.colentries[xCellFaceOrientations.colstart[cell]+j-1] = 1
                else
                    xCellFaceOrientations[j, cell] = 1
                end
            else
                for k = 1 : nfacenodes
                    facenodes[nfacenodes + 1 - k] = xFaceNodes[k,face]
                end
                found_configuration = false
                n = 0
                while !found_configuration
                    n += 1
                    if facenodes[n] == xCellNodes[face_rule[j,1],cell]
                        found_configuration = true
                    end
                end
                n = mod(n-1,3) + 1
                if singleEG == false
                    xCellFaceOrientations.colentries[xCellFaceOrientations.colstart[cell]+j-1] = 1+n
                else
                    xCellFaceOrientations[j, cell] = 1+n
                end
            end
        end

    end
    return xCellFaceOrientations
end

# FaceEdges = Edges for each face (in 3D)
function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{FaceEdges}) where {Tc,Ti}
    xFaceNodes = xgrid[FaceNodes]
    xFaceCells = xgrid[FaceCells]
    xEdgeNodes = xgrid[EdgeNodes]
    xCellEdges = xgrid[CellEdges]
    xFaceEdges = VariableTargetAdjacency(Ti)
    xFaceGeometries = xgrid[FaceGeometries]
    nfaces::Ti = num_sources(xFaceNodes)
    nnodes::Ti = size(xgrid[Coordinates],2)

    # find unique edge enumeration rules
    EG = unique(xFaceGeometries)
    edge_rules = Array{Array{Ti,2},1}(undef,length(EG))
    maxedgenodes = 0
    for j = 1 : length(EG)
        edge_rules[j] = face_enum_rule(EG[j]) # face edges are faces of the face
        maxedgenodes = max(size(edge_rules[j],2),maxedgenodes)
    end
    edge_rule::Array{Ti,2} = edge_rules[1]
    if length(EG) == 1
        singleEG = true
        xFaceEdges = zeros(Ti,nfaces_for_geometry(EG[1]),nfaces)
    else
        singleEG = false
        faceEG = xCellGeometries[1]
        for face = 1 : nfaces
            faceEG = xCellGeometries[cell]
            append!(xCellEdges,zeros(Ti,nfaces_for_geometry(cellEG)))
        end   
    end

    nfacenodes::Ti = num_targets(xFaceNodes,1)
    ncelledges::Ti = 0
    nedgenodes::Ti = 0
    nfaceedges::Ti = nedges_for_geometry(EG[1])
    node::Ti = 0
    cell::Ti = 0
    edge::Ti = 0
    faceedge::Ti = 0
    edge_is_in_face::Bool = false
    found_pos::Bool = false
    flag4face = zeros(Bool,nnodes)
    flag4edge = zeros(Bool,nnodes)
    faceEG = EG[1]
    iEG::Ti = 1
    pos::Ti = 0
    for face = 1 : nfaces

        if singleEG == false
            faceEG = xFaceGeometries[face]
            nfaceedges = nedges_for_geometry(faceEG)
            # find EG index for geometry
            for j=1:length(EG)
                if faceEG == EG[j]
                    iEG = j
                    break;
                end
            end
            edge_rule = edge_rules[iEG] # determines local enumeration of face edges
            nfacenodes = num_targets(xFaceNodes,face)
        end

        # mark nodes of face
        for j = 1 : nfacenodes
            node = xFaceNodes[j, face]
            flag4face[node] = true; 
        end


        # find edges in first adjacent cell
        cell = xFaceCells[1,face]
        ncelledges = num_targets(xCellEdges,cell)
        faceedge = 0
        for cedge = 1 : ncelledges
            edge = xCellEdges[cedge,cell]
            nedgenodes = num_targets(xEdgeNodes,edge)
            edge_is_in_face = true
            for j = 1 : nedgenodes
                if flag4face[xEdgeNodes[j, edge]] == false
                    edge_is_in_face = false
                    break;
                end
            end

            if edge_is_in_face

                # register edge
                # it is important obey same local ordering as in edge_rule
                # to ensure correct dof handling (e.g. for P2-FEM on boundary)

                # mark nodes of edge
                for j = 1 : nedgenodes
                    node = xEdgeNodes[j, edge]
                    flag4edge[node] = true; 
                end

                pos = 0
                found_pos = false
                while found_pos == false
                    pos += 1
                    found_pos = true
                    for k = 1 : nedgenodes
                        if flag4edge[xFaceNodes[edge_rule[pos,k], face]] == false
                            found_pos = false
                            break;
                        end
                    end
                end
                if singleEG
                    xFaceEdges[pos,face] = edge
                else
                    xFaceEdges.colentries[xFaceEdges.colstart[face]+pos-1] = edge
                end

                # reset flag4edge
                for j = 1 : nedgenodes
                    node = xEdgeNodes[j, edge]
                    flag4edge[node] = false; 
                end

                faceedge += 1
                if faceedge == nfaceedges
                    break;
                end
            end
        end

        #reset flag4face
        for j = 1 : nfacenodes
            node = xFaceNodes[j, face]
            flag4face[node] = false; 
        end

    end
    
    xFaceEdges
end

# CellFaces = faces for each cell
function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{CellFaces}) where {Tc,Ti}
    ExtendableGrids.instantiate(xgrid, FaceNodes)
    xgrid[CellFaces]
end

function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{FaceGeometries}) where {Tc,Ti}
    ExtendableGrids.instantiate(xgrid, FaceNodes)
    xgrid[FaceGeometries]
end

# CellEdges = edges for each cell
function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{CellEdges}) where {Tc,Ti}
    ExtendableGrids.instantiate(xgrid, EdgeNodes)
    xgrid[CellEdges]
end

# CellEdgeSigns = edges for each cell
function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{CellEdgeSigns}) where {Tc,Ti}
    ExtendableGrids.instantiate(xgrid, EdgeNodes)
    xgrid[CellEdgeSigns]
end

# CellFaceSigns = orientation signs for each face on each cell
function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{CellFaceSigns}) where {Tc,Ti}
    ExtendableGrids.instantiate(xgrid, FaceNodes)
    xgrid[CellFaceSigns]
end


function collectVolumes4Geometries(T::Type{<:Real}, xgrid::ExtendableGrid{Tc,Ti}, ItemType) where {Tc,Ti}
    # get links to other stuff
    xCoordinates = xgrid[Coordinates]
    xCoordinateSystem = xgrid[CoordinateSystem]
    xItemNodes = xgrid[GridComponent4TypeProperty(ItemType,PROPERTY_NODES) ]
    xGeometries = xgrid[GridComponent4TypeProperty(ItemType,PROPERTY_GEOMETRY) ]
    EG = xgrid[GridComponent4TypeProperty(ItemType,PROPERTY_UNIQUEGEOMETRY) ]
    nitems::Ti = num_sources(xItemNodes)

    # get Volume4ElemType handlers
    handlers = Array{Function,1}(undef, length(EG))
    for j = 1: length(EG)
        handlers[j] = Volume4ElemType(xCoordinates, xItemNodes, EG[j], xCoordinateSystem)
    end

    # init Volumes
    xVolumes::Array{T,1} = zeros(T,nitems)

    # loop over items and call handlers
    iEG::Ti = 1
    itemEG = EG[1]
    for item = 1 : nitems
        if length(EG) > 1
            itemEG = xGeometries[item]
            for j=1:length(EG)
                if itemEG == EG[j]
                    iEG = j
                    break;
                end
            end
        end
        xVolumes[item] = handlers[iEG](item)
    end

    xVolumes
end


function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{CellVolumes}) where {Tc,Ti}
    collectVolumes4Geometries(Float64, xgrid, ITEMTYPE_CELL)
end

function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{FaceVolumes}) where {Tc,Ti}
    collectVolumes4Geometries(Float64, xgrid, ITEMTYPE_FACE)
end

function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{BFaceVolumes}) where {Tc,Ti}
    collectVolumes4Geometries(Float64, xgrid, ITEMTYPE_BFACE)
end

function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{EdgeVolumes}) where {Tc,Ti}
    collectVolumes4Geometries(Float64, xgrid, ITEMTYPE_EDGE)
end

function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{BEdgeVolumes}) where {Tc,Ti}
    collectVolumes4Geometries(Float64, xgrid, ITEMTYPE_BEDGE)
end


function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{BFaces}) where {Tc,Ti}
    # get links to other stuff
    xCoordinates = xgrid[Coordinates]
    xFaceNodes = xgrid[FaceNodes]
    xBFaceNodes = xgrid[BFaceNodes]
    nnodes::Ti = num_sources(xCoordinates)
    nbfaces::Ti = num_sources(xBFaceNodes)

    # init BFaces
    xBFaces::Array{Ti,1} = zeros(Ti,nbfaces)
    #xBFaceGeometries = xgrid[BFaceGeometries]
    #if typeof(xBFaceGeometries) == VectorOfConstants{DataType}
    #    EG = xBFaceGeometries[1]
    #    xBFaceGeometries = Array{DataType,1}(undef,nbfaces)
    #    for j = 1 : nbfaces
    #        xBFaceGeometries[j] = EG
    #    end
    #end

    # transpose FaceNodes to get NodeFaces
    xNodeFaces = atranspose(xFaceNodes)

    flag4item::Array{Bool,1} = zeros(Bool,nnodes)
    nodes_per_bface::Ti = 0
    nodes_per_face::Ti = 0
    common_nodes::Ti = 0
    node::Ti = 0
    nneighbours::Ti = 0
    for bface = 1 : nbfaces
        nodes_per_bface = num_targets(xBFaceNodes,bface)
        for j = 1 : nodes_per_bface
            node = xBFaceNodes[j,bface]
            flag4item[node] = true
        end    

        # get faces for last node of bface
        nneighbours = num_targets(xNodeFaces,node)

        # loop over faces and find the one that matches the bface
        for n = 1 : nneighbours
            face = xNodeFaces[n,node]
            nodes_per_face = num_targets(xFaceNodes,face)
            common_nodes = 0
            for k = 1 : nodes_per_face
                if flag4item[xFaceNodes[k,face]] == true
                    common_nodes += 1
                else
                    break  
                end
            end          
            if common_nodes == nodes_per_face
                xBFaces[bface] = face
                break
            end
        end

        if xBFaces[bface] == 0
            println("WARNING(BFaces): found no matching face for bface $bface with nodes $(xBFaceNodes[:,bface])")
        end

        for j = 1 : nodes_per_bface
            flag4item[xBFaceNodes[j,bface]] = false
        end    
    end

    # enforce that BFaceNodes have same ordering as FaceNodes
    newBFaceNodes = deepcopy(xBFaceNodes)
    for bface = 1 : nbfaces
        nodes_per_face = num_targets(xBFaceNodes,bface)
        for j = 1 : nodes_per_bface
            if typeof(xBFaceNodes) <: VariableTargetAdjacency
                newBFaceNodes.colentries[newBFaceNodes.colstart[bface]+j-1] = xFaceNodes[j,xBFaces[bface]]
            else
                newBFaceNodes[j,bface] = xFaceNodes[j,xBFaces[bface]]
            end
        end
    end
    xgrid[BFaceNodes] = newBFaceNodes

   # xgrid[BFaceGeometries] = xBFaceGeometries
    xBFaces
end



function xinstantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{BEdgeNodes}) where {Tc,Ti}

    dim = size(xgrid[Coordinates],1)
    if dim==1
        xgrid[BEdges]=zeros(Ti,0)
        xgrid[BEdgeNodes]=zeros(Ti,0,0)
        return xgrid[BEdgeNodes]
    end
    
    if dim==2
        xgrid[BEdges] = xgrid[BFaces]
        xgrid[BEdgeNodes] = xgrid[BFaceNodes]
        return xgrid[BEdgeNodes]
    end

    xBFaces = xgrid[BFaces]
    xEdgeNodes = xgrid[EdgeNodes]
    xFaceEdges = xgrid[FaceEdges]
    xBFaceGeometries = xgrid[BFaceGeometries]
    nbfaces = length(xBFaces)
    xBEdges = zeros(Ti,0)

    EG = Triangle2D
    edge::Ti = 0
    face::Ti = 0
    nfaceedges::Ti = 0

    for bface = 1 : nbfaces
        EG = xBFaceGeometries[bface]
        nfaceedges = nedges_for_geometry(EG)
        face = xBFaces[bface]
        for k = 1 : nfaceedges
            edge = xFaceEdges[k,face]
            if !(edge in xBEdges)
                push!(xBEdges, edge)
            end

        end
    end

    nbedges = length(xBEdges)
    xBEdgeNodes = zeros(Ti,2,nbedges)
    for bedge = 1 : nbedges, k = 1 : 2
        xBEdgeNodes[k,bedge] = xEdgeNodes[k,xBEdges[bedge]]
    end

    xgrid[BEdges] = xBEdges
    xBEdgeNodes
end


function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{BEdges}) where {Tc,Ti}
    ExtendableGrids.instantiate(xgrid, BEdgeNodes)
    xgrid[BEdges]
end


function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{FaceCells}) where {Tc,Ti}
    ExtendableGrids.instantiate(xgrid, FaceNodes)
    xgrid[FaceCells]
end

function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{FaceRegions}) where {Tc,Ti}
    # interior faces get region number 0, boundary faces get their boundary region
    xBFaces = xgrid[BFaces]
    xBFaceRegions = xgrid[BFaceRegions]
    xFaceRegions = zeros(Ti,num_sources(xgrid[FaceNodes]))
    for j = 1 : length(xBFaces)
        xFaceRegions[xBFaces[j]] = xBFaceRegions[j]
    end
    xFaceRegions
end



### Not a good idea to have them as vectors of constants
# function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{EdgeRegions}) where {Tc,Ti}
#     return VectorOfConstants(Ti(0),num_sources(xgrid[EdgeNodes]))
# end

# function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{BEdgeRegions}) where {Tc,Ti}
#     return VectorOfConstants(Ti(0),num_sources(xgrid[BEdgeNodes]))
# end

# function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{BEdgeGeometries}) where {Tc,Ti}
#     return VectorOfConstants(Edge1D,num_sources(xgrid[BEdgeNodes]))
# end



function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{BFaceCellPos}) where {Tc,Ti}
    # gets the local position of a bface in the CellFaces row of its adjacent cell

    # get links to other stuff
    xCoordinates = xgrid[Coordinates]
    xCellFaces = xgrid[CellFaces]
    xFaceCells = xgrid[FaceCells]
    xBFaces = xgrid[BFaces]
    nbfaces = length(xBFaces)

    # init BFaces
    xBFaceCellPos = zeros(Ti,nbfaces)

    cface = 0
    cell = 0
    nfaces4cell = 0
    for bface = 1 : nbfaces
        cface = xBFaces[bface]
        cell = xFaceCells[1,cface]
        nfaces4cell = num_targets(xCellFaces,cell)
        for face = 1 : nfaces4cell
            if cface == xCellFaces[face,cell]
                xBFaceCellPos[bface] = face
                break
            end
        end
    end

    xBFaceCellPos
end


function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{FaceNormals}) where {Tc,Ti}
    dim = size(xgrid[Coordinates],1) 
    xCoordinates = xgrid[Coordinates]
    xFaceNodes = xgrid[FaceNodes]
    nfaces = num_sources(xFaceNodes)
    xFaceGeometries = xgrid[FaceGeometries]
    xCoordinateSystem = xgrid[CoordinateSystem]
    xFaceNormals::Array{Float64,2} = zeros(Float64,dim,nfaces)
    normal = zeros(Float64,dim)
    for face = 1 : nfaces
        Normal4ElemType!(normal,xCoordinates,xFaceNodes,face,xFaceGeometries[face],xCoordinateSystem)
        for k = 1 : dim
            xFaceNormals[k, face] = normal[k]
        end    
    end

    xFaceNormals
end


function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{EdgeTangents}) where {Tc,Ti}
    dim = size(xgrid[Coordinates],1) 
    xCoordinates = xgrid[Coordinates]
    xEdgeNodes = xgrid[EdgeNodes]
    nedges = num_sources(xEdgeNodes)
    xEdgeGeometries = xgrid[EdgeGeometries]
    xCoordinateSystem = xgrid[CoordinateSystem]
    xEdgeTangents = zeros(Float64,dim,nedges)
    tangent = zeros(Float64,dim)
    for edge = 1 : nedges
        Tangent4ElemType!(tangent,xCoordinates,xEdgeNodes,edge,xEdgeGeometries[edge],xCoordinateSystem)
        for k = 1 : dim
            xEdgeTangents[k, edge] = tangent[k]
        end    
    end

    xEdgeTangents
end

function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{UniqueCellGeometries}) where {Tc,Ti}
    xUniqueCellGeometries = unique(xgrid[CellGeometries])
end

function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{UniqueFaceGeometries}) where {Tc,Ti}
    xUniqueFaceGeometries = unique(xgrid[FaceGeometries])
end

function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{UniqueBFaceGeometries}) where {Tc,Ti}
    xUniqueBFaceGeometries = unique(xgrid[BFaceGeometries])
end

function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{UniqueEdgeGeometries}) where {Tc,Ti}
    xUniqueEdgeGeometries = unique(xgrid[EdgeGeometries])
end

function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tc,Ti}, ::Type{UniqueBEdgeGeometries}) where {Tc,Ti}
    xUniqueBEdgeGeometries = unique(xgrid[BEdgeGeometries])
end
