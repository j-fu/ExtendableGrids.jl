module ExtendableGrids

using Triangulate
using DocStringExtensions
using ElasticArrays
using AbstractTrees
using InteractiveUtils
using Printf

include("adjacency.jl")
export Adjacency,VariableTargetAdjacency,FixedTargetAdjacency
export atranspose,num_targets,num_sources,num_links,append!, max_num_targets_per_source



include("vectorofconstants.jl")
export VectorOfConstants

include("typehierarchy.jl")
export AbstractExtendableGridApexType
export typehierarchy

include("elementgeometry.jl")
export AbstractElementGeometry, ElementInfo
export elementgeometries

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
export coordinatesystems
export AbstractCoordinateSystem
export Cartesian1D,Cartesian2D,Cartesian3D
export Cylindrical2D,Cylindrical3D
export Polar2D,Polar1D ,Spherical3D,Spherical1D  



include("extendablegrid.jl")
export ExtendableGrid
export instantiate, veryform
export AbstractGridComponent
export AbstractGridAdjacency,AbstractElementGeometries,AbstractElementRegions
export Coordinates,CellNodes,BFaceNodes,CellGeometries,BFaceGeometries,CellRegions,BFaceRegions
export NumCellRegions,NumBFaceRegions,CoordinateSystem
export AbstractGridFloatArray1D,AbstractGridFloatArray2D
export AbstractGridIntegerArray1D,AbstractGridIntegerArray2D
export index_type, coord_type
export dim_space, dim_grid
export num_nodes, num_cells, num_bfaces
export gridcomponents
 

include("subgrid.jl")
export subgrid

include("regionedit.jl")
export cellmask!,bfacemask!

include("simplexgrid.jl")
export simplexgrid, geomspace,glue
export XCoordinates, YCoordinates

include("gridfactory.jl")
export GridFactory
export point!
export facet!
export cellregion!
export flags!

include("tokenstream.jl")
export TokenStream, gettoken, expecttoken,trytoken

include("plot.jl")
export plot
export isplots,isvtkview,ispyplot
export tridata,rectdata

end # module
