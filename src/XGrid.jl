module XGrid
# -> ExtendableGrids


using Triangulate
using DocStringExtensions

include("adjacency.jl")
export Adjacency,VariableTargetAdjacency,FixedTargetAdjacency
export atranspose,num_targets,num_sources,num_links,append!, max_num_targets_per_source

include("vectorofconstants.jl")
export VectorOfConstants

include("elementgeometry.jl")
export AbstractElementGeometry, ElementInfo

export AbstractElementGeometry0D
export Vertex0D

export AbstractElementGeometry1D
export Edge1D

export AbstractElementGeometry2D
export Polygon2D,Triangle2D,Quadrilateral2D,Pentagon2D,Hexagon2D,Parallelogram2D,Circle2D

export AbstractElementGeometry3D
export Polyhedron3D,Tetrahedron3D, Hexahedron3D,Parallelepiped3D,Prism3D,TrianglePrism3D,Sphere3D

export AbstractElementGeometry4D
export Polychoron4D,HyperCube4D

export dim_element


include("coordinatesystem.jl")

export Cartesian1D,Cartesian2D,Cartesian3D
export Cylindrical2D,Cylindrical3D
export Polar2D,Polar1D ,Spherical3D,Spherical1D  



include("extendablegrid.jl")
export ExtendableGrid
export instantiate, veryform
export AbstractGridComponent
export AbstractGridAdjacency,AbstractElementTypes,AbstractElementRegions
export Coordinates,CellNodes,BFaceNodes,CellTypes,BFaceTypes,CellRegions,BFaceRegions
export NumCellRegions,NumBFaceRegions,CoordinateSystem
export AbstractGridFloatArray1D,AbstractGridFloatArray2D
export AbstractGridIntegerArray1D,AbstractGridIntegerArray2D
export index_type, coord_type
export dim_space, dim_grid
export num_nodes, num_cells, num_bfaces

    
include("subgrid.jl")
export subgrid

include("regionedit.jl")
export cellmask!,facemask!

include("simplexgrid.jl")
export simplexgrid, geomspace,glue

include("plot.jl")
export plot

end # module
