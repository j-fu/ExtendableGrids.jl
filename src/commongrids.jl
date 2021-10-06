#####################
# REFERENCE DOMAINS #
#####################

"""
````
    reference_domain(EG::Type{<:AbstractElementGeometry}, T::Type{<:Real} = Float64; scale = [1,1,1], shift = [0,0,0]) -> ExtendableGrid{T,Int32}
````

Generates an ExtendableGrid{T,Int32} for the reference domain of the specified Element Geometry. With scale and shift the coordinates can be manipulated.
"""
function reference_domain(EG::Type{<:AbstractElementGeometry}, T::Type{<:Real} = Float64; scale = [1,1,1], shift = [0,0,0])
    xgrid=ExtendableGrid{T,Int32}()
    xCoordinates=Array{T,2}(refcoords_for_geometry(EG))
    for j = 1 : size(xCoordinates,1)
        xCoordinates[j,:] .+= shift[j]
        xCoordinates[j,:] .*= scale[j]
    end
    xgrid[Coordinates] = xCoordinates
    xCellNodes = zeros(Int32,num_nodes(EG),1)
    xCellNodes[:] = 1:num_nodes(EG)
    xgrid[CellNodes] = xCellNodes
    xgrid[CellGeometries] = VectorOfConstants(EG,1);
    xgrid[CellRegions]=ones(Int32,1)
    xgrid[BFaceRegions]=Array{Int32,1}(1:num_faces(EG))
    xgrid[BFaceNodes]=Array{Int32,2}(local_cellfacenodes(EG))
    xgrid[BFaceGeometries]=VectorOfConstants(facetype_of_cellface(EG, 1), num_faces(EG))
    if dim_element(EG) == 0
        xgrid[CoordinateSystem]=Cartesian0D
    elseif dim_element(EG) == 1
        xgrid[CoordinateSystem]=Cartesian1D
    elseif dim_element(EG) == 2
        xgrid[CoordinateSystem]=Cartesian2D
    elseif dim_element(EG) == 3
        xgrid[CoordinateSystem]=Cartesian3D
    end
    return xgrid
end

"""
````
    Triangle(coords) -> ExtendableGrid{T,Int32}
````

Generates a single triangle with the given coordinates, that should be a 2 x 3 array
with the coordinates of the three vertices, e.g. coords = [0.0 0.0; 1.0 0.0; 0.0 1.0]'.
"""
function grid_triangle(coords::AbstractArray{T,2}) where {T}
    xgrid = reference_domain(Triangle2D)
    xgrid[Coordinates]=coords
    return xgrid
end

# unit cube as one cell with six boundary regions (bottom, front, right, back, left, top)
function grid_unitcube(EG::Type{<:Hexahedron3D}; scale = [1,1,1], shift = [0,0,0])
    return reference_domain(EG; scale = scale, shift = shift)
end

# unit cube as six tets with six boundary regions (bottom, front, right, back, left, top)
function grid_unitcube(::Type{Tetrahedron3D}; scale = [1,1,1], shift = [0,0,0])
    xgrid=ExtendableGrid{Float64,Int32}()
    xCoordinates=Array{Float64,2}([0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1]')
    for j = 1 : size(xCoordinates,1)
        xCoordinates[j,:] .+= shift[j]
        xCoordinates[j,:] .*= scale[j]
    end
    xgrid[Coordinates] = xCoordinates

    xCellNodes=Array{Int32,2}([1 2 3 7; 1 3 4 7; 1 5 6 7; 1 8 5 7; 1 6 2 7;1 4 8 7]')
    xgrid[CellNodes] = xCellNodes
    xgrid[CellGeometries] = VectorOfConstants(Tetrahedron3D,6);
    ncells = num_sources(xCellNodes)
    xgrid[CellRegions]=ones(Int32,ncells)
    xgrid[BFaceRegions]=Array{Int32,1}([1,1,2,2,3,3,4,4,5,5,6,6])
    xBFaceNodes=Array{Int32,2}([1 3 2; 1 4 3; 1 2 6;1 6 5;2 3 7;2 7 6;3 4 7;7 4 8;8 4 1;1 5 8; 5 6 7; 5 7 8]')
    xgrid[BFaceNodes]=xBFaceNodes
    nbfaces = num_sources(xBFaceNodes)
    xgrid[BFaceGeometries]=VectorOfConstants(Triangle2D,nbfaces)
    xgrid[CoordinateSystem]=Cartesian3D
    return xgrid
end



# unit square as one cell with four boundary regions (bottom, right, top, left)
function grid_unitsquare(EG::Type{<:Quadrilateral2D}; scale = [1,1], shift = [0,0])
    return reference_domain(EG; scale = scale, shift = shift)
end

# unit square as two triangles with four boundary regions (bottom, right, top, left)
function grid_unitsquare(::Type{<:Triangle2D}; scale = [1,1], shift = [0,0])
    xgrid=ExtendableGrid{Float64,Int32}()
    xCoordinates=Array{Float64,2}([0 0; 1 0; 1 1; 0 1; 0.5 0.5]')
    for j = 1 : size(xCoordinates,1)
        xCoordinates[j,:] .+= shift[j]
        xCoordinates[j,:] .*= scale[j]
    end
    xgrid[Coordinates] = xCoordinates
    xgrid[CellNodes]=Array{Int32,2}([1 2 5; 2 3 5; 3 4 5; 4 1 5]')
    xgrid[CellGeometries]=VectorOfConstants(Triangle2D,4)
    xgrid[CellRegions]=ones(Int32,4)
    xgrid[BFaceRegions]=Array{Int32,1}([1,2,3,4])
    xgrid[BFaceNodes]=Array{Int32,2}([1 2; 2 3; 3 4; 4 1]')
    xgrid[BFaceGeometries]=VectorOfConstants(Edge1D,4)
    xgrid[CoordinateSystem]=Cartesian2D
    return xgrid
end


# unit square as two triangles with four boundary regions (bottom, right, top, left)
function grid_lshape(::Type{<:Triangle2D}; scale = [1,1], shift = [0,0])
    xgrid=ExtendableGrid{Float64,Int32}()
    xCoordinates=Array{Float64,2}([0 0; 1 0; 1 1; 0 1; -1 1; -1 0; -1 -1; 0 -1]')
    for j = 1 : size(xCoordinates,1)
        xCoordinates[j,:] .+= shift[j]
        xCoordinates[j,:] .*= scale[j]
    end
    xgrid[Coordinates] = xCoordinates
    xgrid[CellNodes]=Array{Int32,2}([4 2 3; 2 4 1; 6 4 5; 4 6 1; 6 8 1; 8 6 7]')
    xgrid[CellGeometries]=VectorOfConstants(Triangle2D,6)
    xgrid[CellRegions]=Array{Int32,1}([1,1,1,1,1,1])
    xgrid[BFaceRegions]=Array{Int32,1}([1,2,3,4,5,6,7,8])
    xgrid[BFaceNodes]=Array{Int32,2}([1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 1]')
    xgrid[BFaceGeometries]=VectorOfConstants(Edge1D,8)
    xgrid[CoordinateSystem]=Cartesian2D
    return xgrid
end

# unit suqare as mixed triangles and squares with four boundary regions (bottom, right, top, left)
function grid_unitsquare_mixedgeometries()

    xgrid=ExtendableGrid{Float64,Int32}()
    xgrid[Coordinates]=Array{Float64,2}([0 0; 4//10 0; 1 0; 0 6//10; 4//10 6//10; 1 6//10;0 1; 4//10 1; 1 1]')
    xCellNodes=VariableTargetAdjacency(Int32)
    xCellGeometries=[Triangle2D, Triangle2D, Parallelogram2D, Parallelogram2D, Triangle2D, Triangle2D];
    
    append!(xCellNodes,[1,5,4])
    append!(xCellNodes,[1,2,5])
    append!(xCellNodes,[2,3,6,5])
    append!(xCellNodes,[4,5,8,7]) 
    append!(xCellNodes,[5,6,9])
    append!(xCellNodes,[8,5,9])

    xgrid[CellNodes] = xCellNodes
    xgrid[CellGeometries] = xCellGeometries
    ncells = num_sources(xCellNodes)
    xgrid[CellRegions]=ones(Int32,ncells)
    xgrid[BFaceRegions]=Array{Int32,1}([1,1,2,2,3,3,4,4])
    xBFaceNodes=Array{Int32,2}([1 2; 2 3; 3 6; 6 9; 9 8; 8 7; 7 4; 4 1]')
    xgrid[BFaceNodes]=xBFaceNodes
    nbfaces = num_sources(xBFaceNodes)
    xgrid[BFaceGeometries]=VectorOfConstants(Edge1D,nbfaces)
    xgrid[CoordinateSystem]=Cartesian2D

    return xgrid
end

