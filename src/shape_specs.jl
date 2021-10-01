#########################################
# VERTEX0D <: AbstractElementGeometry0D #       [1]
#########################################

refcoords_for_geometry(::Type{<:AbstractElementGeometry0D}) = [0]
nnodes_for_geometry(::Type{<:AbstractElementGeometry0D}) = 1
nfaces_for_geometry(::Type{<:AbstractElementGeometry0D}) = 0
nedges_for_geometry(::Type{<:AbstractElementGeometry0D}) = 0


#######################################      
# EDGE1D <: AbstractElementGeometry1D #     [1]-----[2]        [1] = 0, [2] = 1
#######################################      

refcoords_for_geometry(::Type{<:AbstractElementGeometry1D}) = [0; 1]

nnodes_for_geometry(::Type{<:AbstractElementGeometry1D}) = 2
nfaces_for_geometry(::Type{<:AbstractElementGeometry1D}) = 2
nedges_for_geometry(::Type{<:AbstractElementGeometry1D}) = 0
face_enum_rule(::Type{<:AbstractElementGeometry1D}) = reshape([1; 2],2,1)
facetype_of_cellface(P1::Type{<:AbstractElementGeometry1D},P2::Type{<:AbstractElementGeometry1D}, k) = Vertex0D
facetype_of_cellface(::Type{<:AbstractElementGeometry1D}, k) = Vertex0D

xrefFACE2xrefCELL(::Type{<:AbstractElementGeometry1D}) = [ [(xref4FACE) -> [1]],
                                                           [(xref4FACE) -> [1]] ]

xrefFACE2xrefOFACE(::Type{<:AbstractElementGeometry1D}) = [(xref4FACE) -> xref4FACE, (xref4FACE) -> 1 .- xref4FACE]

#                   [3]                 
#                    | \   
##############       |   \                    [1] = (0,0)
# Triangle2D #       |     \                  [2] = (1,0)
##############       |       \                [3] = (0,1)
#                    |         \ 
#                   [1]--------[2]

refcoords_for_geometry(::Type{<:Triangle2D}) = [0 0; 1 0; 0 1]

nnodes_for_geometry(::Type{<:Triangle2D}) = 3
face_enum_rule(::Type{<:Triangle2D}) = [1 2; 2 3; 3 1]
edge_enum_rule(T::Type{<:Triangle2D}) = face_enum_rule(T)

# maps of reference coords on cell face to reference coords in cell
xrefFACE2xrefCELL(::Type{<:Triangle2D}) = [ (xref4FACE) -> [xref4FACE[1],0],
                                            (xref4FACE) -> [1-xref4FACE[1],xref4FACE[1]],
                                            (xref4FACE) -> [0,1-xref4FACE[1]]
                                            ]

# maps of reference coords on face to reference coords in face with other orientation
xrefFACE2xrefOFACE(::Type{<:Triangle2D}) = [(xref4FACE) -> xref4FACE,                                    # orientation 1 = [1,2,3]
                                            (xref4FACE) -> [xref4FACE[1],1-xref4FACE[1]-xref4FACE[2]],   # orientation 2 = [1,3,2]
                                            (xref4FACE) -> [1-xref4FACE[1]-xref4FACE[2],xref4FACE[2]],   # orientation 3 = [3,2,1]
                                            (xref4FACE) -> [xref4FACE[2],xref4FACE[1]]                   # orientation 4 = [2,1,3]
                                            ]   

#                        [4]--------[3]               
#                         |          |             [1] = (0,0)
###################       |          |             [2] = (1,0)
# Quadrilateral2D #       |          |             [3] = (1,1)
###################       |          |             [4] = (0,1)
#                        [1]--------[2]

refcoords_for_geometry(::Type{<:Quadrilateral2D}) = [0 0; 1 0; 1 1; 0 1]

nnodes_for_geometry(::Type{<:Quadrilateral2D}) = 4
face_enum_rule(::Type{<:Quadrilateral2D}) = [1 2; 2 3; 3 4; 4 1]
edge_enum_rule(T::Type{<:Quadrilateral2D}) = face_enum_rule(T)

# maps of reference coords on cell face to reference coords in cell
xrefFACE2xrefCELL(::Type{<:Quadrilateral2D}) = [ (xref4FACE) -> [xref4FACE[1],0],
                                                 (xref4FACE) -> [1,xref4FACE[1]],
                                                 (xref4FACE) -> [1-xref4FACE[1],1],
                                                 (xref4FACE) -> [0,1-xref4FACE[1]]
                                                 ]


#############################      
# AbstractElementGeometry2D #    
#############################   

facetype_of_cellface(P1::Type{<:AbstractElementGeometry2D}, k) = Edge1D
facetype_of_cellface(P1::Type{<:AbstractElementGeometry2D},P2::Type{<:AbstractElementGeometry2D}, k) = Edge1D
nfaces_for_geometry(EG::Type{<:AbstractElementGeometry2D}) = nnodes_for_geometry(EG)
nedges_for_geometry(EG::Type{<:AbstractElementGeometry2D}) = nnodes_for_geometry(EG)



#                      [4]                 
#                       |\\   
#################       | \ \                    [1] = (0,0,0)
# Tetrahedron3D #       |  \  \                  [2] = (1,0,0)
#################       |   \   \                [3] = (0,1,0)
#                       | _-[3]-_ \              [4] = (0,0,1)
#                      [1]--------[2]

refcoords_for_geometry(::Type{<:Tetrahedron3D}) = [0 0 0; 1 0 0; 0 1 0; 0 0 1]

nfaces_for_geometry(::Type{<:Tetrahedron3D}) = 4
nnodes_for_geometry(::Type{<:Tetrahedron3D}) = 4
nedges_for_geometry(::Type{<:Tetrahedron3D}) = 6
face_enum_rule(::Type{<:Tetrahedron3D}) = [1 3 2; 1 2 4; 2 3 4; 1 4 3]
facetype_of_cellface(P1::Type{<:Tetrahedron3D},P2::Type{<:Tetrahedron3D}, k) = Triangle2D
facetype_of_cellface(::Type{<:Tetrahedron3D}, k) = Triangle2D
edge_enum_rule(::Type{<:Tetrahedron3D}) = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]
celledges_for_cellface(::Type{<:Tetrahedron3D}) = [2 4 1; 1 5 3; 4 6 5; 3 6 2]

# maps of reference coords on cell face to reference coords in cell
xrefFACE2xrefCELL(::Type{<:Tetrahedron3D}) = [ (xref4FACE) -> [xref4FACE[2],xref4FACE[1],0],
                                               (xref4FACE) -> [xref4FACE[1],0,xref4FACE[2]],
                                               (xref4FACE) -> [1-xref4FACE[1]-xref4FACE[2],xref4FACE[1],xref4FACE[2]],
                                               (xref4FACE) -> [0,xref4FACE[2],xref4FACE[1]]                 
                                                ]


#                         
#                         [8]--------[7]
#                        / |        / |          [1] = (0,0,0)
#                     [5]--------[6]  |          [2] = (1,0,0)
################       |   |      |   |          [3] = (1,1,0)
# Hexahedron3D #       |   |      |   |          [4] = (0,1,0)
################       |   |      |   |          [5] = (0,0,1)
#                      |  [4]-----|--[3]         [6] = (1,0,1)
#                      | /        | /            [7] = (1,1,1)
#                     [1]--------[2]             [8] = (0,1,1)

refcoords_for_geometry(::Type{<:Hexahedron3D}) = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1]

nnodes_for_geometry(::Type{<:Hexahedron3D}) = 8
nfaces_for_geometry(::Type{<:Hexahedron3D}) = 6
nedges_for_geometry(::Type{<:Hexahedron3D}) = 12
face_enum_rule(::Type{<:Hexahedron3D}) = [4 3 2 1; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 5 6 7 8]
facetype_of_cellface(P1::Type{<:Hexahedron3D},P2::Type{<:Hexahedron3D}, k) = Quadrilateral2D
facetype_of_cellface(P1::Type{<:Parallelepiped3D},P2::Type{<:Hexahedron3D}, k) = Parallelogram2D
facetype_of_cellface(::Type{<:Hexahedron3D}, k) = Quadrilateral2D
facetype_of_cellface(::Type{<:Parallelepiped3D}, k) = Parallelogram2D
edge_enum_rule(::Type{<:Hexahedron3D}) = [1 2; 2 3; 3 4; 4 1; 1 5; 2 6; 3 7; 4 8; 5 6; 6 7; 7 8; 8 5]
celledges_for_cellface(::Type{<:Hexahedron3D}) = [3 2 1 4; 1 6 9 5; 2 7 10 6; 3 8 11 7; 4 5 12 8; 9 10 11 12]


#############################      
# AbstractElementGeometry3D #    
#############################  

edgetype_of_celledge(::Type{<:AbstractElementGeometry2D}, k) = Edge1D
edgetype_of_celledge(::Type{<:AbstractElementGeometry3D}, k) = Edge1D




###############
### VOLUMES ###
###############

function Volume4ElemType(Coords, Nodes, ::Type{<:Vertex0D}, ::Type{<:ExtendableGrids.AbstractCoordinateSystem})
    function closure(item)
        return 0
    end
end

function Volume4ElemType(Coords, Nodes, ::Type{<:Edge1D}, ::Type{Cartesian1D})
    function closure(item)
        return abs(Coords[1, Nodes[2,item]] - Coords[1, Nodes[1,item]])
    end
end

function Volume4ElemType(Coords, Nodes, ::Type{<:Edge1D}, ::Type{Cartesian2D})
    function closure(item)
        return sqrt((Coords[1, Nodes[2,item]] - Coords[1, Nodes[1,item]]).^2 + (Coords[2, Nodes[2,item]] - Coords[2, Nodes[1,item]]).^2)
    end
end

function Volume4ElemType(Coords, Nodes, ::Type{<:Edge1D}, ::Type{Cartesian3D})
    function closure(item)
        return sqrt((Coords[1, Nodes[2,item]] - Coords[1, Nodes[1,item]]).^2 + (Coords[2, Nodes[2,item]] - Coords[2, Nodes[1,item]]).^2  + (Coords[3, Nodes[2,item]] - Coords[3, Nodes[1,item]]).^2)
    end
end

function Volume4ElemType(Coords, Nodes, ::Type{<:Triangle2D}, ::Type{Cartesian2D})
    function closure(item)
        return 1 // 2 * ( Coords[1, Nodes[1, item]] * (Coords[2, Nodes[2,item]] -  Coords[2, Nodes[3, item]])
                      +   Coords[1, Nodes[2, item]] * (Coords[2, Nodes[3,item]] -  Coords[2, Nodes[1, item]])
                      +   Coords[1, Nodes[3, item]] * (Coords[2, Nodes[1,item]] -  Coords[2, Nodes[2, item]]) )
    end
end

function Volume4ElemType(Coords, Nodes, ::Type{<:Quadrilateral2D}, ::Type{Cartesian2D})
    function closure(item)
        return 1//2 * (   (Coords[1, Nodes[1, item]] - Coords[1, Nodes[3, item]]) * (Coords[2, Nodes[2, item]] - Coords[2, Nodes[4, item]])
                    + (Coords[1, Nodes[4, item]] - Coords[1, Nodes[2, item]]) * (Coords[2, Nodes[1, item]] - Coords[2, Nodes[3, item]]) );

    end
end

function Volume4ElemType(Coords, Nodes, ::Type{<:Triangle2D}, ::Type{Cartesian3D})
    d12 = zeros(Float64,3)
    d14 = zeros(Float64,3)
    function closure(item)
        # norm(cross(p(1)-p(2), p(1)-p(3)), 2)
        for k = 1 : 3
            d12[k] = Coords[k, Nodes[1, item]] - Coords[k, Nodes[2, item]]
            d14[k] = Coords[k, Nodes[1, item]] - Coords[k, Nodes[3, item]]
        end
        return sqrt((d12[2]*d14[3]-d12[3]*d14[2])^2 + (d12[3]*d14[1]-d12[1]*d14[3])^2 + (d12[1]*d14[2]-d12[2]*d14[1])^2) / 2;
    end
end


function Volume4ElemType(Coords, Nodes, ::Type{<:Parallelogram2D}, ::Type{Cartesian3D})
    d12 = zeros(Float64,3)
    d14 = zeros(Float64,3)
    function closure(item)
        # norm(cross(p(1)-p(2), p(1)-p(3)), 2)
        for k = 1 : 3
            d12[k] = Coords[k, Nodes[1, item]] - Coords[k, Nodes[2, item]]
            d14[k] = Coords[k, Nodes[1, item]] - Coords[k, Nodes[4, item]]
        end
        return sqrt((d12[2]*d14[3]-d12[3]*d14[2])^2 + (d12[3]*d14[1]-d12[1]*d14[3])^2 + (d12[1]*d14[2]-d12[2]*d14[1])^2);
    end
end


function Volume4ElemType(Coords, Nodes, ::Type{<:Parallelepiped3D}, ::Type{Cartesian3D})
    function closure(item)
        return    ((Coords[1, Nodes[5, item]] - Coords[1, Nodes[1, item]]) * ( (Coords[2, Nodes[2, item]] - Coords[2, Nodes[1, item]]) * (Coords[3, Nodes[4, item]] - Coords[3, Nodes[1, item]]) - (Coords[2, Nodes[4, item]] - Coords[2, Nodes[1, item]]) * (Coords[3, Nodes[2, item]] - Coords[3, Nodes[1, item]])) 
        + (Coords[2, Nodes[5, item]] - Coords[2, Nodes[1, item]]) * ( (Coords[3, Nodes[2, item]] - Coords[3, Nodes[1, item]]) * (Coords[1, Nodes[4, item]] - Coords[1, Nodes[1, item]]) - (Coords[1, Nodes[2, item]] - Coords[1, Nodes[1, item]]) * (Coords[3, Nodes[4, item]] - Coords[3, Nodes[1, item]])) 
        + (Coords[3, Nodes[5, item]] - Coords[3, Nodes[1, item]]) * ( (Coords[1, Nodes[2, item]] - Coords[1, Nodes[1, item]]) * (Coords[2, Nodes[4, item]] - Coords[2, Nodes[1, item]]) - (Coords[2, Nodes[2, item]] - Coords[2, Nodes[1, item]]) * (Coords[1, Nodes[4, item]] - Coords[1, Nodes[1, item]])));
    
    end
end


function Volume4ElemType(Coords, Nodes, ::Type{<:Tetrahedron3D}, ::Type{Cartesian3D})
    function closure(item)
        return    1 // 6 * ((Coords[1, Nodes[4, item]] - Coords[1, Nodes[1, item]]) * ( (Coords[2, Nodes[2, item]] - Coords[2, Nodes[1, item]]) * (Coords[3, Nodes[3, item]] - Coords[3, Nodes[1, item]]) - (Coords[2, Nodes[3, item]] - Coords[2, Nodes[1, item]]) * (Coords[3, Nodes[2, item]] - Coords[3, Nodes[1, item]])) 
        + (Coords[2, Nodes[4, item]] - Coords[2, Nodes[1, item]]) * ( (Coords[3, Nodes[2, item]] - Coords[3, Nodes[1, item]]) * (Coords[1, Nodes[3, item]] - Coords[1, Nodes[1, item]]) - (Coords[1, Nodes[2, item]] - Coords[1, Nodes[1, item]]) * (Coords[3, Nodes[3, item]] - Coords[3, Nodes[1, item]])) 
        + (Coords[3, Nodes[4, item]] - Coords[3, Nodes[1, item]]) * ( (Coords[1, Nodes[2, item]] - Coords[1, Nodes[1, item]]) * (Coords[2, Nodes[3, item]] - Coords[2, Nodes[1, item]]) - (Coords[2, Nodes[2, item]] - Coords[2, Nodes[1, item]]) * (Coords[1, Nodes[3, item]] - Coords[1, Nodes[1, item]])));
    end
end
  

###############
### NORMALS ###
###############

function Normal4ElemType!(normal, Coords, Nodes, item, ::Type{<:Vertex0D}, ::Type{Cartesian2D})
    normal[1] = 0.0
    normal[2] = 0.0
end

function Normal4ElemType!(normal, Coords, Nodes, item, ::Type{<:Vertex0D}, ::Type{Cartesian1D})
    normal[1] = 1.0
end

function Normal4ElemType!(normal, Coords, Nodes, item, ::Type{<:Edge1D}, ::Type{Cartesian2D})
    # rotate tangent
    normal[1] = Coords[2, Nodes[2,item]] - Coords[2, Nodes[1,item]]
    normal[2] = Coords[1, Nodes[1,item]] - Coords[1, Nodes[2,item]]
    # divide by length
    normal ./= sqrt(normal[1]^2+normal[2]^2)
end

function Normal4ElemType!(normal, Coords, Nodes, item, ::Type{<:Quadrilateral2D}, ::Type{Cartesian3D})
    # cross(p(1)-p(2), p(1)-p(4)) / length
    normal[1]  = (Coords[2, Nodes[1, item]] - Coords[2, Nodes[2, item]]) * (Coords[3, Nodes[1, item]] - Coords[3, Nodes[4, item]])
    normal[1] -= (Coords[3, Nodes[1, item]] - Coords[3, Nodes[2, item]]) * (Coords[2, Nodes[1, item]] - Coords[2, Nodes[4, item]])
    normal[2]  = (Coords[3, Nodes[1, item]] - Coords[3, Nodes[2, item]]) * (Coords[1, Nodes[1, item]] - Coords[1, Nodes[4, item]])
    normal[2] -= (Coords[1, Nodes[1, item]] - Coords[1, Nodes[2, item]]) * (Coords[3, Nodes[1, item]] - Coords[3, Nodes[4, item]])
    normal[3]  = (Coords[1, Nodes[1, item]] - Coords[1, Nodes[2, item]]) * (Coords[2, Nodes[1, item]] - Coords[2, Nodes[4, item]])
    normal[3] -= (Coords[2, Nodes[1, item]] - Coords[2, Nodes[2, item]]) * (Coords[1, Nodes[1, item]] - Coords[1, Nodes[4, item]])
    # divide by length
    normal ./= sqrt(normal[1]^2+normal[2]^2+normal[3]^2)

    ## old version
    #d12 = @views Coords[:, Nodes[1, item]] - Coords[:, Nodes[2, item]]
    #d14 = @views Coords[:, Nodes[1, item]] - Coords[:, Nodes[4, item]]
    #normal[1] = d12[2]*d14[3]-d12[3]*d14[2]
    #normal[2] = d12[3]*d14[1]-d12[1]*d14[3]
    #normal[3] = d12[1]*d14[2]-d12[2]*d14[1]
end

function Normal4ElemType!(normal, Coords, Nodes, item, ::Type{<:Triangle2D}, ::Type{Cartesian3D})
    # cross(p(1)-p(2), p(1)-p(3)) / length
    normal[1]  = (Coords[2, Nodes[1, item]] - Coords[2, Nodes[2, item]]) * (Coords[3, Nodes[1, item]] - Coords[3, Nodes[3, item]])
    normal[1] -= (Coords[3, Nodes[1, item]] - Coords[3, Nodes[2, item]]) * (Coords[2, Nodes[1, item]] - Coords[2, Nodes[3, item]])
    normal[2]  = (Coords[3, Nodes[1, item]] - Coords[3, Nodes[2, item]]) * (Coords[1, Nodes[1, item]] - Coords[1, Nodes[3, item]])
    normal[2] -= (Coords[1, Nodes[1, item]] - Coords[1, Nodes[2, item]]) * (Coords[3, Nodes[1, item]] - Coords[3, Nodes[3, item]])
    normal[3]  = (Coords[1, Nodes[1, item]] - Coords[1, Nodes[2, item]]) * (Coords[2, Nodes[1, item]] - Coords[2, Nodes[3, item]])
    normal[3] -= (Coords[2, Nodes[1, item]] - Coords[2, Nodes[2, item]]) * (Coords[1, Nodes[1, item]] - Coords[1, Nodes[3, item]])
    # divide by length
    normal ./= sqrt(normal[1]^2+normal[2]^2+normal[3]^2)
end


################
### TANGENTS ###
################

function Tangent4ElemType!(tangent, Coords, Nodes, item, ::Type{<:Edge1D}, ::Type{Cartesian2D})
    tangent[1] = Coords[1,Nodes[2,item]] - Coords[1, Nodes[1,item]]
    tangent[2] = Coords[2,Nodes[2,item]] - Coords[2, Nodes[1,item]]
    # divide by length
    tangent ./= sqrt(tangent[1]^2+tangent[2]^2)
end

function Tangent4ElemType!(tangent, Coords, Nodes, item, ::Type{<:Edge1D}, ::Type{Cartesian3D})
    tangent[1] = Coords[1,Nodes[2,item]] - Coords[1, Nodes[1,item]]
    tangent[2] = Coords[2,Nodes[2,item]] - Coords[2, Nodes[1,item]]
    tangent[3] = Coords[3,Nodes[2,item]] - Coords[3, Nodes[1,item]]
    # divide by length
    tangent ./= sqrt(tangent[1]^2+tangent[2]^2+tangent[3]^2)
end
