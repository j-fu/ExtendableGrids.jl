module XGrid
# -> ExtendableGrids
using Triangulate

include("adjacency.jl")
export atranspose,ntargets,nsources,nlinks,append!
export Adjacency,VariableTargetAdjacency,FixedTargetAdjacency

include("vectorofconstants.jl")
export VectorOfConstants

include("elementinfo.jl")
export AbstractElementType, ElementInfo
export Simplex1D, Simplex2D, Simplex3D

include("extendablegrid.jl")
export ExtendableGrid
export instantiate, veryform
export AbstractGridComponent
export AbstractGridArray1D,AbstractGridArray2D, AbstractGridAdjacency,AbstractElementTypes,AbstractElementRegions
export Coordinates,CellNodes,BFaceNodes,CellTypes,BFaceTypes,CellRegions,BFaceRegions
    
include("generate.jl")
export generate

include("plot.jl")
export plot

end # module
