function bulk_mark(xgrid, refinement_indicators, theta = 0.5; indicator_AT = ON_CELLS)

    xFaceCells = xgrid[FaceCells]
    xFaceNodes = xgrid[FaceNodes]
    nfaces = num_sources(xFaceNodes)

    # redistribute refinement_indicators to face indicators (if necessary)
    if indicator_AT == ON_CELLS
        refind4face = zeros(Float64,nfaces)
        cell::Int = 0
        for face = 1 : nfaces
            for k = 1 : 2
                cell = xFaceCells[k,face]
                if cell > 0
                    refind4face[face] += refinement_indicators[cell]
                end
            end
            if cell > 0
                refind4face[face] /= 2
            end
        end
    elseif indicator_AT == ON_FACES
        refind4face = refinement_indicators
    end

    # mark largest indicators until theta*total is reached
    p = Base.sortperm(refind4face[:], rev = true)
    totalsum = sum(refind4face)
    csum = 0
    j = 0
    facemarker = zeros(Bool,nfaces)
    while csum <= theta*totalsum
        j += 1
        csum += refind4face[p[j]]
        facemarker[p[j]] = true
    end

    return facemarker
end



# adaptive refinement rules fo TRIANGLE2D
# first k nodes are the CellNodes
# next m nodes are the CellFaces midpoints
# next node is the CellMidpoint (if needed)
red_refinement(::Type{<:Triangle2D}) = [1 4 6; 4 2 5; 6 5 3; 5 6 4]
blueR_refinement(::Type{<:Triangle2D}) = [3 1 4; 4 2 5; 3 4 5]
blueL_refinement(::Type{<:Triangle2D}) = [1 4 6; 2 3 4; 4 3 6]
green_refinement(::Type{<:Triangle2D}) = [3 1 4; 2 3 4]

function RGB_refinement_rule(EG::Type{<:Triangle2D}, marked_faces::Array{Bool,1})
    if any(marked_faces) == true
        @assert marked_faces[1] = true # closure should always mark reference face
        if marked_faces[2] == true && marked_faces[3] == true
            return red_refinement(EG), 1
        elseif marked_faces[2] == true && marked_faces[3] == false
            return blueR_refinement(EG), 2
        elseif marked_faces[2] == false && marked_faces[3] == true
            return blueL_refinement(EG), 3
        elseif marked_faces[2] == false && marked_faces[3] == false
            return green_refinement(EG), 4
        end
    else
        return Array{Int,2}([1 2 3]), 5
    end
end



"""
$(TYPEDSIGNATURES)

generates a new ExtendableGrid by red-green-blue mesh refinement of triangular meshes, see e.g.

Carstensen, C.
--An Adaptive Mesh-Refining Algorithm Allowing for an H^1 Stable L^2 Projection onto Courant Finite Element Spaces--
Constr Approx 20, 549–564 (2004). https://doi.org/10.1007/s00365-003-0550-5

The bool array facemarkers determines which faces should be bisected. Note, that a closuring is performed
such that the first face in every triangle with a marked face is also refined.
"""
function RGB_refine(source_grid::ExtendableGrid{T,K}, facemarkers::Array{Bool,1}) where {T,K}
    
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
    
    xCellNodes = VariableTargetAdjacency(Int32)
    xCellGeometries = []
    xCellRegions = zeros(Int32,0)
    oldCellNodes = source_grid[CellNodes]
    oldCellFaces = source_grid[CellFaces]
    oldFaceNodes = source_grid[FaceNodes]
    nfaces = num_sources(oldFaceNodes)
    ncells = num_sources(oldCellNodes)
    nrefcounts = [0,0,0,0,0]

    # closuring
    # @logmsg MoreInfo "RGB refinement with $(sum(facemarkers)) marked faces"
    is_refined = true
    closure_finished = false
    while closure_finished == false
        closure_finished = true
        for cell = 1 : ncells
            is_refined = false
            for f = 1 : 3
                if facemarkers[oldCellFaces[f,cell]] == true
                    is_refined = true
                    break;
                end
            end
            if is_refined && facemarkers[oldCellFaces[1,cell]] == false
                # mark also reference edge
                facemarkers[oldCellFaces[1,cell]] = true
                closure_finished = false
            end
        end
    end
    # @logmsg DeepInfo "marked faces after closure = $(sum(facemarkers))"


    # determine number of new vertices
    oldvertices = size(oldCoordinates,2)
    newvertices = 0
    for face = 1 : nfaces
        if facemarkers[face] == true
            newvertices += 1
        end
    end
    xCoordinates = zeros(Float64,size(oldCoordinates,1),oldvertices+newvertices)
    @views xCoordinates[:,1:oldvertices] = oldCoordinates

    # refine faces
    newvertex = zeros(Float64,size(xCoordinates,1))
    newnodenr = oldvertices
    newnode4facemidpoint = zeros(Int,nfaces)
    nnodes4item = 0
    # add face midpoints to Coordinates
    for face = 1 : nfaces
        if facemarkers[face] == true
            nnodes4item = num_targets(oldFaceNodes,face)
            fill!(newvertex,0.0)
            for k = 1 : nnodes4item, d = 1 : size(xCoordinates,1)
                newvertex[d] += xCoordinates[d,oldFaceNodes[k,face]] 
            end    
            newvertex ./= nnodes4item
            newnodenr += 1
            newnode4facemidpoint[face] = newnodenr
            for d = 1 : size(xCoordinates,1)
                xCoordinates[d,newnodenr] = newvertex[d]
            end
        end
    end    
    
    # determine new cells
    nnodes4item = 3
    nfaces4item = 3
    ncells = 0
    iEG = 1
    subitemnodes = zeros(Int32,max_num_targets_per_source(oldCellNodes)+max_num_targets_per_source(oldCellFaces)+1)
    m = 0
    localmarked_faces = zeros(Bool,3)
    refine_rule = Array{Int,2}([1 2 3])
    for cell = 1 : num_sources(oldCellNodes)
        for f = 1 : 3
            localmarked_faces[f] = facemarkers[oldCellFaces[f,cell]]
        end
        refine_rule, refid = RGB_refinement_rule(Triangle2D,localmarked_faces)
        nrefcounts[refid] += 1
        for j = 1 : size(refine_rule,1)
            for k = 1 : size(refine_rule,2)
                m = refine_rule[j,k]
                if m <= nnodes4item 
                    subitemnodes[k] = oldCellNodes[m,cell]
                elseif m <= nnodes4item + nfaces4item
                    subitemnodes[k] = newnode4facemidpoint[oldCellFaces[m-nnodes4item,cell]]
                end       
            end    
            append!(xCellNodes,subitemnodes[1:size(refine_rule,2)])
            push!(xCellGeometries,Triangle2D)
            push!(xCellRegions,oldCellRegions[cell])
        end    
        ncells += size(refine_rule,1)
    end

    # @logmsg DeepInfo "\tred/blueR/blueL/green/unrefined = $nrefcounts"

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
    
    xBFaceRegions = zeros(Int32,0)
    xBFaceGeometries = []
    nbfaces = num_sources(oldBFaceNodes)
    xBFaceNodes = []

    EG = Base.unique(oldBFaceGeometries)

    refine_rules = Array{Array{Int,2},1}(undef,length(EG))
    for j = 1 : length(EG)
        refine_rules[j] = uniform_refine_rule(EG[j])
    end

    newbfaces = 0
    for bface = 1 : nbfaces
        face = oldBFaces[bface]
        itemEG = oldBFaceGeometries[bface]
        if facemarkers[face] == true
            nnodes4item = num_nodes(itemEG)
            nfaces4item = num_faces(itemEG)
            iEG = findfirst(isequal(itemEG), EG)

            for j = 1 : size(refine_rules[iEG],1)
                for k = 1 : size(refine_rules[iEG],2)
                    m = refine_rules[iEG][j,k]
                    if dim == 2
                        if m <= nnodes4item 
                            subitemnodes[k] = oldBFaceNodes[m,bface]
                        else
                            subitemnodes[k] = newnode4facemidpoint[face]
                        end            
                    end
                end
                append!(xBFaceNodes,subitemnodes[1:size(refine_rules[iEG],2)])
                push!(xBFaceGeometries,itemEG)
                push!(xBFaceRegions,oldBFaceRegions[bface])
            end    
            newbfaces +=1
        else
            append!(xBFaceNodes,oldBFaceNodes[:,bface])
            push!(xBFaceGeometries,itemEG)
            push!(xBFaceRegions,oldBFaceRegions[bface])
        end
    end
    @debug "bisected bfaces = $newbfaces"
    
    xgrid[BFaceNodes] = Array{Int32,2}(reshape(xBFaceNodes,(2,nbfaces+newbfaces)))
    xgrid[BFaceRegions] = xBFaceRegions
    xgrid[BFaceGeometries] = Array{DataType,1}(xBFaceGeometries)

    xgrid = split_grid_into(xgrid,Triangle2D)
    return xgrid
end
