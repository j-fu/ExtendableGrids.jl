module XGrid
# -> ExtendableGrids
using Triangulate
using DocStringExtensions

include("adjacency.jl")
export atranspose,num_targets,num_sources,num_links,append!
export Adjacency,VariableTargetAdjacency,FixedTargetAdjacency

include("vectorofconstants.jl")
export VectorOfConstants

include("elementinfo.jl")
export AbstractElementType, ElementInfo
export Simplex0D, Simplex1D, Simplex2D, Simplex3D
export Cartesian1D,Cartesian2D,Cartesian3D
export Cylindrical2D,Cylindrical3D
export Polar2D,Polar1D ,Spherical3D,Spherical1D  

 
include("extendablegrid.jl")
export ExtendableGrid
export instantiate, veryform
export AbstractGridComponent
export AbstractGridArray1D,AbstractGridArray2D, AbstractGridAdjacency,AbstractElementTypes,AbstractElementRegions
export Coordinates,CellNodes,BFaceNodes,CellTypes,BFaceTypes,CellRegions,BFaceRegions
export NumCellRegions,NumBFaceRegions,CoordinateSystem

include("generate.jl")
export generate, simplexgrid

include("plot.jl")
export plot


end # module
