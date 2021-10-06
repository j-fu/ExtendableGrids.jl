

# refinements can store parent informations
abstract type CellParents <: AbstractGridIntegerArray1D end

function ExtendableGrids.instantiate(xgrid::ExtendableGrid, ::Type{CellParents})
    ncells = num_sources(xgrid[CellNodes]) 
    VectorOfConstants{Int32}(0,ncells)
end



# functions that tell how to split one ElementGeometry into another
split_rule(::Type{Triangle2D}, ::Type{Triangle2D}) = reshape([1,2,3],1,3)
split_rule(::Type{Tetrahedron3D}, ::Type{Tetrahedron3D}) = reshape([1,2,3,4],1,4)
split_rule(::Type{<:Quadrilateral2D}, ::Type{Triangle2D}) = [1 2 3;1 3 4]
split_rule(::Type{Edge1D}, ::Type{Triangle2D}) = reshape([1,2,2],1,3)
split_rule(::Type{<:Hexahedron3D}, ::Type{Tetrahedron3D}) = [1 2 4 8; 2 3 4 8; 2 7 3 8; 6 8 7 2; 6 5 8 2; 8 5 1 2]
#[1 2 3 7; 1 3 4 7; 1 5 6 7; 1 8 5 7; 1 6 2 7; 1 4 8 7]

"""
$(TYPEDSIGNATURES)

generates a new ExtendableGrid by splitting each cell into subcells of the specified targetgeometry

split rules exist for
- Quadrilateral2D into Triangle2D
- Hexahedron3D into Tetrahedron3D
"""
function split_grid_into(source_grid::ExtendableGrid{T,K}, targetgeometry::Type{<:AbstractElementGeometry}) where {T,K}
    xgrid=ExtendableGrid{T,K}()
    xgrid[Coordinates]=source_grid[Coordinates]
    oldCellGeometries = source_grid[CellGeometries]
    oldCellRegions = source_grid[CellRegions]
    EG = Base.unique(oldCellGeometries)
    
    split_rules = Array{Array{Int,2},1}(undef,length(EG))
    for j = 1 : length(EG)
        split_rules[j] = split_rule(EG[j],targetgeometry)
    end
    xCellNodes=[]
    xCellRegions = zeros(Int32,0)
    oldCellNodes=source_grid[CellNodes]
    nnodes4item = 0
    ncells = 0
    itemEG = targetgeometry
    iEG = 1
    for cell = 1 : num_sources(oldCellNodes)
        nnodes4item = num_targets(oldCellNodes,cell)
        itemEG = oldCellGeometries[cell]
        iEG = findfirst(isequal(itemEG), EG)
        for j = 1 : size(split_rules[iEG],1), k = 1 : size(split_rules[iEG],2)
            append!(xCellNodes,oldCellNodes[split_rules[iEG][j,k],cell])
        end    
        for j = 1 : size(split_rules[iEG],1)
            push!(xCellRegions,oldCellRegions[cell])
        end
        ncells += size(split_rules[iEG],1)
    end
    xCellNodes = reshape(xCellNodes,num_nodes(targetgeometry),ncells)
    xgrid[CellNodes] = Array{Int32,2}(xCellNodes)
    xgrid[CellGeometries] = VectorOfConstants(targetgeometry,ncells)
    if typeof(oldCellRegions) <: VectorOfConstants
        xgrid[CellRegions] = VectorOfConstants{Int32}(1,ncells)
    else
        xgrid[CellRegions] = xCellRegions
    end

    # find new boundary faces (easy in 2D, not so easy in 3D)
    if dim_element(targetgeometry) == 2 # BFaces are Edge1D wich stay the same
        xgrid[BFaceNodes]=source_grid[BFaceNodes]
        xgrid[BFaceRegions]=source_grid[BFaceRegions]
        xgrid[BFaceGeometries]=VectorOfConstants(facetype_of_cellface(targetgeometry,1),num_sources(xgrid[BFaceNodes]))
    elseif dim_element(targetgeometry) == 3 
        # BFaces may be split into different shapes, e.g. from Quadrilateral2D to two Triangle2D
        # and it is hard to predict how they are splitted
        # so we do something lazy here and search for new faces that lie in old bfaces
        oldBFaceNodes = source_grid[BFaceNodes]
        oldBFaceRegions = source_grid[BFaceRegions]
        newFaceNodes = xgrid[FaceNodes]
        nfaces = num_sources(newFaceNodes) 
        nbfaces = num_sources(oldBFaceNodes)
        
        newBFaceNodes = []
        newBFaceRegions = []
        newnbfaces = 0
        nnodes = size(xgrid[Coordinates],2)
        flag4item = zeros(Bool,nnodes)
        nodes_per_bface::Int = 0
        nodes_per_face::Int = 0
        common_nodes::Int = 0
        for bface = 1 : nbfaces
            # flag nodes of old bface
            nodes_per_bface = num_targets(oldBFaceNodes,bface)
            for j = 1 : nodes_per_bface
                flag4item[oldBFaceNodes[j,bface]] = true
            end    
            
            # find matching faces
            for face = 1 : nfaces
                nodes_per_face = num_targets(newFaceNodes,face)
                common_nodes = 0
                for k = 1 : nodes_per_face
                    if flag4item[newFaceNodes[k,face]] == true
                        common_nodes += 1
                    else
                        break  
                    end
                end          
                if common_nodes == nodes_per_face
                    append!(newBFaceNodes,newFaceNodes[:,face])
                    push!(newBFaceRegions,oldBFaceRegions[bface])
                    newnbfaces += 1
                end
            end

            # reset flags
            for j = 1 : nodes_per_bface
                flag4item[oldBFaceNodes[j,bface]] = false
            end    
        end

        newBFaceNodes = reshape(newBFaceNodes,num_nodes(facetype_of_cellface(targetgeometry,1)),newnbfaces)
        xgrid[BFaceNodes]=Array{Int32,2}(newBFaceNodes)
        xgrid[BFaceRegions]=Array{Int32,1}(newBFaceRegions)
        xgrid[BFaceGeometries]=VectorOfConstants(facetype_of_cellface(targetgeometry,1),newnbfaces)
    end
    xgrid[CoordinateSystem]=source_grid[CoordinateSystem]
    return xgrid
end


uniform_refine_needfacemidpoints(::Type{<:AbstractElementGeometry}) = true
uniform_refine_needcellmidpoints(::Type{<:AbstractElementGeometry}) = false

# uniform refinement rules in 1D
# first k nodes are the CellNodes
# next node is the CellMidpoint
uniform_refine_rule(::Type{<:Edge1D}) = [1 3; 3 2]
uniform_refine_needcellmidpoints(::Type{<:Edge1D}) = true

# uniform refinement rules in 2D
# first k nodes are the CellNodes
# next m nodes are the CellFaces midpoints
# next node is the CellMidpoint (if needed)
uniform_refine_rule(::Type{<:Triangle2D}) = [1 4 6; 4 2 5; 6 5 3; 5 6 4]
uniform_refine_rule(::Type{<:Quadrilateral2D}) = [1 5 9 8; 2 6 9 5; 3 7 9 6; 4 8 9 7]
uniform_refine_needcellmidpoints(::Type{<:Quadrilateral2D}) = true

# uniform refinement rules in 3D
# first k nodes are the CellNodes
# next m nodes are the CellEdges midpoints
# next n nodes are the CellFaces midpoints
# next node is the CellMidpoint (if needed)
uniform_refine_rule(::Type{<:Tetrahedron3D}) = [ 1 5 6 7;
                                                 2 5 9 8;
                                                 3 10 6 8;
                                                 4 10 9 7;
                                                 10 5 8 9;
                                                 5 10 7 9;
                                                 5 10 8 6;
                                                 10 5 7 6]
uniform_refine_needfacemidpoints(::Type{<:Tetrahedron3D}) = false
uniform_refine_rule(::Type{<:Hexahedron3D}) = [     1   9  21  12  13  22  27  25;
                                                    9   2  10  21  22  14  23  27;
                                                   12  21  11   4  25  27  24  16;
                                                   21  10   3  11  27  23  15  24;
                                                   13  22  27  25   5  17  26  20;
                                                   22  14  23  27  17   6  18  26;
                                                   25  27  24  16  20  26  19   8;
                                                   27  23  15  24  26  18   7  19]
uniform_refine_needcellmidpoints(::Type{<:Hexahedron3D}) = true


"""
$(TYPEDSIGNATURES)

generates a new ExtendableGrid by uniform refinement of each cell in the given grid

uniform refinement rules are available for these AbstractElementGeometries:
- Line1D (bisection into two subsegments)
- Triangle2D (red refinement into four subtriangles)
- Quadrilateral2D (into four subquadrilaterals)
- Tetrahedron (into eight subtetrahedrons)
- Hexahedron (into eight subhexahedrons)

if multiple geometries are in the mesh uniform refinement will only work
if all refinement rules refine faces and edges (in 3D) equally
(so no hanging nodes are created)
"""
function uniform_refine(source_grid::ExtendableGrid{T,K}; store_parents = false) where {T,K}
    # @logmsg MoreInfo "Uniform refinement of $(num_sources(source_grid[CellNodes])) cells"
    
    xgrid = ExtendableGrid{T,K}()
    xgrid[CoordinateSystem]=source_grid[CoordinateSystem]

    # unpack stuff from source grid
    oldCoordinates = source_grid[Coordinates]
    oldCellGeometries = source_grid[CellGeometries]
    oldCellRegions = source_grid[CellRegions]
    EG = Base.unique(oldCellGeometries)

    # get dimension of CellGeometries
    # currently it is assumed to be the same for all cells
    dim = dim_element(EG[1]) 
    
    refine_rules = Array{Array{Int,2},1}(undef,length(EG))
    need_facemidpoints = uniform_refine_needfacemidpoints(EG[1])
    for j = 1 : length(EG)
        refine_rules[j] = uniform_refine_rule(EG[j])
        @assert uniform_refine_needfacemidpoints(EG[j]) == need_facemidpoints
    end

    xCellNodes = VariableTargetAdjacency(Int32)
    xCellGeometries = []
    xCellRegions = zeros(Int32,0)
    oldCellNodes = source_grid[CellNodes]
    oldCellFaces = source_grid[CellFaces]
    oldCellEdges = []
    nfaces = 0
    if dim > 1
        oldFaceNodes = source_grid[FaceNodes]
        nfaces = num_sources(oldFaceNodes)
    end
    nedges = 0
    if dim > 2 
        oldEdgeNodes = source_grid[EdgeNodes]
        nedges = num_sources(oldEdgeNodes)
        oldCellEdges = source_grid[CellEdges]
    end


    # determine number of new vertices
    itemEG = Triangle2D
    newvertices = 0 # in 1D no additional vertices on the faces are needed
    if dim == 2 # in 2D each face is halved
        newvertices = nfaces
    elseif dim == 3 # in 2D each face and edge is halved
        newvertices = nedges
        if need_facemidpoints 
            newvertices += nfaces
        end
    end
    oldvertices = size(oldCoordinates,2)
    newnode = oldvertices + newvertices
    # additionally cell midpoints are needed for some refinements
    for cell = 1 : num_sources(oldCellNodes)
        itemEG = oldCellGeometries[cell]
        if uniform_refine_needcellmidpoints(itemEG) == true
            newvertices += 1
        end    
    end
    xCoordinates = zeros(Float64,size(oldCoordinates,1),oldvertices+newvertices)
    @views xCoordinates[:,1:oldvertices] = oldCoordinates

    
    newvertex = zeros(Float64,size(xCoordinates,1))
    nnodes4item = 0
    if dim > 2 # add edge midpoints to Coordinates
        for edge = 1 : nedges
            nnodes4item = num_targets(oldEdgeNodes,edge)
            fill!(newvertex,0.0)
            for k = 1 : nnodes4item, d = 1 : size(xCoordinates,1)
                newvertex[d] += xCoordinates[d,oldEdgeNodes[k,edge]] 
            end    
            newvertex ./= nnodes4item
            for d = 1 : size(xCoordinates,1)
                xCoordinates[d,oldvertices+edge] = newvertex[d]
            end
        end    
    end
    if dim > 1 && need_facemidpoints 
        # add face midpoints to Coordinates
        for face = 1 : nfaces
            nnodes4item = num_targets(oldFaceNodes,face)
            fill!(newvertex,0.0)
            for k = 1 : nnodes4item, d = 1 : size(xCoordinates,1)
                newvertex[d] += xCoordinates[d,oldFaceNodes[k,face]] 
            end    
            newvertex ./= nnodes4item
            for d = 1 : size(xCoordinates,1)
                xCoordinates[d,oldvertices+nedges+face] = newvertex[d]
            end
        end    
    end
    
    # determine new cells
    nnodes4item = 0
    nfaces4item = 0
    nedges4item = 0
    ncells = 0
    iEG = 1
    subitemnodes = zeros(Int32,max_num_targets_per_source(oldCellNodes)+max_num_targets_per_source(oldCellFaces)+1)
    m = 0
    xCellParents = zeros(Int32,0)
    for cell = 1 : num_sources(oldCellNodes)
        itemEG = oldCellGeometries[cell]
        nnodes4item = num_nodes(itemEG)
        nfaces4item = num_faces(itemEG)
        nedges4item = num_edges(itemEG)
        iEG = findfirst(isequal(itemEG), EG)
        if uniform_refine_needcellmidpoints(itemEG) == true
            # add cell midpoint to Coordinates
            newnode += 1
            fill!(newvertex,0.0)
            for k = 1 : nnodes4item, d = 1 : size(xCoordinates,1)
                newvertex[d] += xCoordinates[d,oldCellNodes[k,cell]] 
            end    
            newvertex ./= nnodes4item
            for d = 1 : size(xCoordinates,1)
                xCoordinates[d,newnode] = newvertex[d]
            end
        end
        for j = 1 : size(refine_rules[iEG],1)
            for k = 1 : size(refine_rules[iEG],2)
                m = refine_rules[iEG][j,k]
                if dim == 1
                    if m <= nnodes4item 
                        subitemnodes[k] = oldCellNodes[m,cell]
                    else
                        subitemnodes[k] = newnode
                    end
                elseif dim == 2
                    if m <= nnodes4item 
                        subitemnodes[k] = oldCellNodes[m,cell]
                    elseif m <= nnodes4item + nfaces4item
                        subitemnodes[k] = oldvertices + oldCellFaces[m-nnodes4item,cell]
                    else
                        subitemnodes[k] = newnode
                    end        
                elseif dim == 3
                    if m <= nnodes4item 
                        subitemnodes[k] = oldCellNodes[m,cell]
                    elseif m <= nnodes4item + nedges4item
                        subitemnodes[k] = oldvertices + oldCellEdges[m-nnodes4item,cell]
                    elseif m <= nnodes4item + nedges4item + nfaces4item
                        subitemnodes[k] = oldvertices + nedges + oldCellFaces[m-nnodes4item-nedges4item,cell]
                    else
                        subitemnodes[k] = newnode
                    end        
                end
            end    
            append!(xCellNodes,subitemnodes[1:size(refine_rules[iEG],2)])
            push!(xCellGeometries,itemEG)
            push!(xCellRegions, oldCellRegions[cell])
            if store_parents
                push!(xCellParents,cell)
            end
        end    
        ncells += size(refine_rules[iEG],1)
    end

    # assign new cells to grid
    xgrid[Coordinates] = xCoordinates
    if typeof(oldCellNodes) == Array{Int32,2}
        nnodes4item = size(oldCellNodes,1)
        xgrid[CellNodes] = reshape(xCellNodes.colentries,nnodes4item,num_sources(xCellNodes))
    else
        xgrid[CellNodes] = xCellNodes
    end
    if typeof(oldCellRegions) <: VectorOfConstants
        xgrid[CellRegions] = VectorOfConstants{Int32}(1,ncells)
    else
        xgrid[CellRegions] = xCellRegions
    end
    xgrid[CellGeometries] = Array{DataType,1}(xCellGeometries)


    # determine new boundary faces
    oldBFaceNodes = source_grid[BFaceNodes]
    oldBFaces = source_grid[BFaces]
    oldBFaceRegions = source_grid[BFaceRegions]
    oldBFaceGeometries = source_grid[BFaceGeometries]
    oldBFacesCellPos = source_grid[BFaceCellPos]
    oldFaceCells = source_grid[FaceCells]
    
    if dim == 1
        xgrid[BFaceNodes] = oldBFaceNodes
        xgrid[BFaceRegions] = oldBFaceRegions
        xgrid[BFaceGeometries] = oldBFaceGeometries
    else
        xBFaceRegions = zeros(Int32,0)
        xBFaceGeometries = []
        nbfaces = num_sources(oldBFaceNodes)
        if dim == 2 || typeof(oldBFaceNodes) == Array{Int32,2}
            xBFaceNodes = []
        else
            xBFaceNodes = VariableTargetAdjacency(Int32)
        end
        if dim == 3
            xCellEdges = source_grid[CellEdges]
            xNodeEdges = atranspose(oldEdgeNodes)
        else
            xCellEdges = []
        end

        EG = Base.unique(oldBFaceGeometries)

        refine_rules = Array{Array{Int,2},1}(undef,length(EG))
        for j = 1 : length(EG)
            refine_rules[j] = uniform_refine_rule(EG[j])
        end


        bcell = 0
        edge = 0
        newnbfaces = 0
        for bface = 1 : nbfaces
            face = oldBFaces[bface]
            itemEG = oldBFaceGeometries[bface]
            nnodes4item = num_nodes(itemEG)
            nfaces4item = num_faces(itemEG)
            iEG = findfirst(isequal(itemEG), EG)
            bface_enum_rule = local_cellfacenodes(itemEG)

            for j = 1 : size(refine_rules[iEG],1)
                for k = 1 : size(refine_rules[iEG],2)
                    m = refine_rules[iEG][j,k]
                    if dim == 2
                        if m <= nnodes4item 
                            subitemnodes[k] = oldBFaceNodes[m,bface]
                        else
                            subitemnodes[k] = oldvertices + face
                        end        
                    elseif dim == 3
                        if m <= nnodes4item 
                            subitemnodes[k] = oldBFaceNodes[m,bface]
                        elseif m <= nnodes4item + nfaces4item
                            edge = m-nnodes4item # local number
                            # find global edge number
                            edge = intersect(xNodeEdges[:,oldBFaceNodes[bface_enum_rule[1,edge],bface]],xNodeEdges[:,oldBFaceNodes[bface_enum_rule[2,edge],bface]])[1]
                            subitemnodes[k] = oldvertices + edge
                        else
                            subitemnodes[k] = oldvertices + nedges + face
                        end        
                    end
                end
                append!(xBFaceNodes,subitemnodes[1:size(refine_rules[iEG],2)])
                push!(xBFaceGeometries,itemEG)
                push!(xBFaceRegions,oldBFaceRegions[bface])
                newnbfaces += 1
            end    
        end
        if dim == 2 || typeof(oldBFaceNodes) == Array{Int32,2}
            xgrid[BFaceNodes] = Array{Int32,2}(reshape(xBFaceNodes,(size(oldBFaceNodes,1),newnbfaces)))
        else
            xgrid[BFaceNodes] = xBFaceNodes
        end
        xgrid[BFaceRegions] = xBFaceRegions
        xgrid[BFaceGeometries] = Array{DataType,1}(xBFaceGeometries)
    end    

    if store_parents
        xgrid[CellParents] = xCellParents
    end


    return xgrid
end


function uniform_refine(source_grid::ExtendableGrid{T,K}, nrefinements::Int; store_parents = false) where {T,K}
    xgrid = source_grid
    parents = 1:num_sources(xgrid[CellNodes])
    for j=1:nrefinements
        xgrid = uniform_refine(xgrid; store_parents = store_parents)
        if store_parents
            parents = parents[xgrid[CellParents]]
        end
    end
    if store_parents
        xgrid[CellParents] = parents
    end
    return xgrid
end


# barycentric refinement rules
# first k nodes are the CellNodes, k+1-th  node is cell midpoint
barycentric_refine_rule(::Type{<:Triangle2D}) = [1 2 4; 2 3 4; 3 1 4]


"""
$(TYPEDSIGNATURES)

generates a new ExtendableGrid by barycentric refinement of each cell in the source grid

barycentric refinement is available for these ElementGeometries
- Quadrilateral2D (first split into Triangle2D)
- Triangle2D
"""
function barycentric_refine(source_grid::ExtendableGrid{T,K}) where {T,K}
    # @logmsg MoreInfo "Barycentric refinement of $(num_sources(source_grid[CellNodes])) cells"
    
    # split first into triangles
    source_grid = split_grid_into(source_grid,Triangle2D)

    xgrid = ExtendableGrid{T,K}()
    oldCoordinates = source_grid[Coordinates]
    oldCellGeometries = source_grid[CellGeometries]
    oldCellRegions = source_grid[CellRegions]
    EG = Base.unique(oldCellGeometries)
    
    refine_rules = Array{Array{Int,2},1}(undef,length(EG))
    for j = 1 : length(EG)
        refine_rules[j] = barycentric_refine_rule(EG[j])
    end
    xCellNodes = VariableTargetAdjacency(Int32)
    xCellGeometries = []
    xCellRegions = zeros(Int32,0)

    oldCellNodes = source_grid[CellNodes]
    oldCellFaces = source_grid[CellFaces]

    # determine number of new vertices
    itemEG = Triangle2D
    newvertices = 0
    for cell = 1 : num_sources(oldCellNodes)
        newvertices += 1
    end
    oldvertices = size(oldCoordinates,2)
    xCoordinates = zeros(Float64,size(oldCoordinates,1),oldvertices+newvertices)
    @views xCoordinates[:,1:oldvertices] = oldCoordinates

    # determine new cells
    nnodes4item = 0
    ncells = 0
    iEG = 1
    subitemnodes = zeros(Int32,max_num_targets_per_source(oldCellNodes)+max_num_targets_per_source(oldCellFaces)+1)
    newnode = oldvertices
    m = 0
    newvertex = zeros(Float64,size(xCoordinates,1))
    for cell = 1 : num_sources(oldCellNodes)
        nnodes4item = num_targets(oldCellNodes,cell)
        nfaces4item = num_targets(oldCellFaces,cell)
        itemEG = oldCellGeometries[cell]
        iEG = findfirst(isequal(itemEG), EG)
        
            # add cell midpoint to Coordinates
            newnode += 1
            fill!(newvertex,0.0)
            for k = 1 : nnodes4item, d = 1 : size(xCoordinates,1)
                newvertex[d] += xCoordinates[d,oldCellNodes[k,cell]] 
            end    
            newvertex ./= nnodes4item
            for d = 1 : size(xCoordinates,1)
                xCoordinates[d,newnode] = newvertex[d]
            end

        for j = 1 : size(refine_rules[iEG],1)
            for k = 1 : size(refine_rules[iEG],2)
                m = refine_rules[iEG][j,k]
                if m <= nnodes4item 
                    subitemnodes[k] = oldCellNodes[m,cell]
                else
                    subitemnodes[k] = newnode
                end        
            end    
            append!(xCellNodes,subitemnodes[1:size(refine_rules[iEG],2)])
            push!(xCellGeometries,itemEG)
            push!(xCellRegions, oldCellRegions[cell])
        end    
        ncells += size(refine_rules[iEG],1)
    end

    xgrid[Coordinates] = xCoordinates
    if typeof(oldCellNodes) == Array{Int32,2}
        nnodes4item = size(oldCellNodes,1)
        xgrid[CellNodes] = reshape(xCellNodes.colentries,nnodes4item,num_sources(xCellNodes))
    else
        xgrid[CellNodes] = xCellNodes
    end
    if typeof(oldCellRegions) <: VectorOfConstants
        xgrid[CellRegions] = VectorOfConstants{Int32}(1,ncells)
    else
        xgrid[CellRegions] = xCellRegions
    end
    xgrid[CellGeometries] = Array{DataType,1}(xCellGeometries)
    xgrid[BFaceNodes]=source_grid[BFaceNodes]
    xgrid[BFaceRegions]=source_grid[BFaceRegions]
    xgrid[BFaceGeometries]=source_grid[BFaceGeometries]
    xgrid[CoordinateSystem]=source_grid[CoordinateSystem]
    return xgrid
end
