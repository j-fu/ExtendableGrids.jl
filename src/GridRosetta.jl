module GridRosetta
using GeometryBasics
using ExtendableGrids

function GeometryBasics.Mesh(grid::ExtendableGrid)
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    npoints=size(coord,2)
    nfaces=size(cellnodes,2)
    points=Vector{Point2f0}(undef,npoints)
    for i=1:npoints
        points[i]=Point2f0(coord[1,i],coord[2,i])
    end
    faces=Vector{GLTriangleFace}(undef,nfaces)
    for i=1:nfaces
        faces[i]=GLTriangleFace(cellnodes[1,i],cellnodes[2,i],cellnodes[3,i])
    end
    Mesh(points,faces)
end


end
