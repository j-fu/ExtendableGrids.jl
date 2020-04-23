module XGrid
# -> ExtendableGrids
using Triangulate
using DocStringExtensions

include("adjacency.jl")
export atranspose,num_targets,num_sources,num_links,append!
export Adjacency,VariableTargetAdjacency,FixedTargetAdjacency

include("vectorofconstants.jl")
export VectorOfConstants

include("elementgeometry.jl")
export AbstractElementGeometry, ElementInfo

export AbstractElementGeometry0D
export AbstractElementGeometry1D
export AbstractElementGeometry2D
export AbstractElementGeometry3D
export AbstractElementGeometry4D

export Vertex0D

export Edge1D

export Polygon2D
export Triangle2D
export Quadrilateral2D
export Pentagon2D
export Hexagon2D
export Parallelogram2D

export Circle2D

export Polyhedron3D
export Tetrahedron3D
export Hexahedron3D
export Parallelepiped3D
export Prism3D
export TrianglePrism3D
export Sphere3D

export Polychoron4D
export HyperCube4D



include("coordinatesystem.jl")

export Cartesian1D,Cartesian2D,Cartesian3D
export Cylindrical2D,Cylindrical3D
export Polar2D,Polar1D ,Spherical3D,Spherical1D  

export dim_element


include("extendablegrid.jl")
export ExtendableGrid
export instantiate, veryform
export AbstractGridComponent
export AbstractGridArray1D,AbstractGridArray2D, AbstractGridAdjacency,AbstractElementTypes,AbstractElementRegions
export Coordinates,CellNodes,BFaceNodes,CellTypes,BFaceTypes,CellRegions,BFaceRegions
export NumCellRegions,NumBFaceRegions,CoordinateSystem
export index_type, coord_type

include("generate.jl")
export generate, simplexgrid

include("subgrid.jl")
export subgrid

include("regionedit.jl")
export cellmask!,facemask!


include("plot.jl")
export plot


end # module
