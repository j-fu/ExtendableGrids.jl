"""
$(TYPEDEF)
Apex type for grid components
"""
abstract type AbstractGridComponent end

"""
$(TYPEDEF)
2D Array on grid components (e.g. coordinates)
"""
abstract type AbstractGridFloatArray2D <: AbstractGridComponent end

"""
$(TYPEDEF)
1D Array on grid components
"""
abstract type AbstractGridFloatArray1D <: AbstractGridComponent end

"""
$(TYPEDEF)
1D Array on grid components
"""
abstract type AbstractGridIntegerArray1D <: AbstractGridComponent end

"""
$(TYPEDEF)
2D Array on grid components
"""
abstract type AbstractGridIntegerArray2D <: AbstractGridComponent end

"""
$(TYPEDEF)
Integer number
"""
abstract type AbstractGridIntegerConstant <: AbstractGridComponent end

"""
$(TYPEDEF)
Floating point  number
"""
abstract type AbstractGridFloatConstant <: AbstractGridComponent end


"""
$(TYPEDEF)
Any kind of adjacency between grid components
"""
abstract type AbstractGridAdjacency <: AbstractGridComponent end


"""
$(TYPEDEF)
ElementInfo{AbstractCellType} on arrays
"""
abstract type AbstractElementTypes <: AbstractGridComponent end

"""
$(TYPEDEF)
ElementInfo{Ti} on arrays
"""
abstract type AbstractElementRegions <: AbstractGridComponent end

"""
$(TYPEDEF)
Basic Grid components with classification via intermediate types
"""
abstract type Coordinates <: AbstractGridFloatArray2D end

"""
$(TYPEDEF)
"""
abstract type CellNodes <: AbstractGridAdjacency end
"""
$(TYPEDEF)
"""
abstract type BFaceNodes <: AbstractGridAdjacency end

"""
$(TYPEDEF)
"""
abstract type CellTypes <: AbstractElementTypes end
"""
$(TYPEDEF)
"""
abstract type BFaceTypes <: AbstractElementTypes end

"""
$(TYPEDEF)
"""
abstract type CellRegions <: AbstractElementRegions end
"""
$(TYPEDEF)
"""
abstract type BFaceRegions <: AbstractElementRegions end

"""
$(TYPEDEF)
"""
abstract type NumCellRegions <: AbstractGridIntegerConstant end
"""
$(TYPEDEF)
"""
abstract type NumBFaceRegions <: AbstractGridIntegerConstant end
"""
$(TYPEDEF)
"""
abstract type CoordinateSystem <: AbstractGridComponent end


"""
$(TYPEDEF)
"""
const ElementInfo{T}=Union{Vector{T},VectorOfConstants{T}}

"""
$(TYPEDEF)
Grid type wrapping Dict
"""
mutable struct ExtendableGrid{Tc,Ti}
    components::Dict{Type{<:AbstractGridComponent},Any}
    ExtendableGrid{Tc,Ti}() where{Tc,Ti} =new(Dict{Type{<:AbstractGridComponent},Any}())
end

"""
$(SIGNATURES)
Default veryform method.

"veryform"  means "verify and/or transform"  and is called to check
and possibly transform components to be added to the grid.


The default method just passes data through.
"""
veryform(grid::ExtendableGrid,v,::Type{<:AbstractGridComponent})=v

"""
$(SIGNATURES)
Check proper type of adjacencies
"""
veryform(grid::ExtendableGrid{Tc,Ti},v,T::Type{<:AbstractGridAdjacency}) where{Tc,Ti}= typeof(v)<:Adjacency{Ti} ? v : throw("Type mismatch")


"""
$(SIGNATURES)
Set new grid component
"""
Base.setindex!(grid::ExtendableGrid,v,T::Type{<:AbstractGridComponent})= grid.components[T]=veryform(grid,v,T)


"""
$(SIGNATURES)
To be called by getindex. This triggers lazy creation of 
non-existing gridcomponents
"""
Base.get!(grid::ExtendableGrid,T::Type{<:AbstractGridComponent})= get!( ()->veryform(grid,instantiate(grid,T),T),  grid.components,T)


"""
$(SIGNATURES)
"Hook" for methods instantiating lazy components. 
See https://white.ucc.asn.au/2020/04/19/Julia-Antipatterns.html
"""
function instantiate end


"""
$(SIGNATURES)
Generic method for obtaining grid component.

This method is mutating in the sense that non-existing grid components
are created on demand.

Due to the fact that components are stored as Any the return
value triggers type instability.
"""
Base.getindex(grid::ExtendableGrid,T::Type{<:AbstractGridComponent})=get!(grid,T)


"""
$(SIGNATURES)
Type specific method to obtain adjacency component
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridAdjacency})::Adjacency{Ti} where{Tc,Ti}
    get!(grid,T)
end

"""
$(SIGNATURES)
Type specific method to obtain 2D array from grid
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridFloatArray2D})::Array{Tc,2} where{Tc,Ti}
    get!(grid,T)
end

"""
$(SIGNATURES)
Type specific method to obtain 1D array from grid
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridFloatArray1D})::Array{Tc,1} where{Tc,Ti}
    get!(grid,T)
end

"""
$(SIGNATURES)
Type specific method to obtain  element type
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractElementTypes})::ElementInfo{DataType} where{Tc,Ti}
    get!(grid,T)
end

"""
$(SIGNATURES)
Type specific method to obtain element region number
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractElementRegions})::ElementInfo{Ti} where{Tc,Ti}
    get!(grid,T)
end


"""
$(SIGNATURES)
Type specific method to obtain  integer constant stored on elements
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridIntegerConstant})::Ti where{Tc,Ti}
    get!(grid,T)
end

"""
$(SIGNATURES)
Type specific method to obtain float constant stored on elements
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridFloatConstant})::Tc where{Tc,Ti}
    get!(grid,T)
end


"""
$(SIGNATURES)
Instantiate number of cell regions
"""
instantiate(grid, ::Type{NumCellRegions})=maximum(grid[CellRegions])
"""
$(SIGNATURES)
Instantiate number of bface regions
"""
instantiate(grid, ::Type{NumBFaceRegions})=maximum(grid[BFaceRegions])


"""
$(SIGNATURES)

Space dimension of grid
"""
dim_space(grid::ExtendableGrid)= size(grid[Coordinates],1)

"""
$(SIGNATURES)

Grid dimension dimension of grid (larges element dimension)
"""
dim_grid(grid::ExtendableGrid)=  dim_element(grid[CellTypes][1])


##########################################################
"""
$(SIGNATURES)

Number of nodes in grid
"""
num_nodes(grid::ExtendableGrid)= size(grid[Coordinates],2)


##########################################################
"""
$(SIGNATURES)

Number of cells in grid
"""
num_cells(grid::ExtendableGrid)= num_sources(grid[CellNodes])

##########################################################
"""
$(SIGNATURES)

Number of edges in grid
"""
num_edges(grid::ExtendableGrid)= num_sources(grid[EdgeNodes])


################################################
"""
$(SIGNATURES)

Number of boundary faces in grid.
"""
num_bfaces(grid::ExtendableGrid)= num_sources(grid[BFaceNodes])

################################################
"""
$(SIGNATURES)

Maximum  cell  region number
"""
num_cellregions(grid::ExtendableGrid)=grid[NumCellRegions]


"""
$(SIGNATURES)

Maximum  boundary face region number
"""
num_bfaceregions(grid::ExtendableGrid)=grid[NumBFaceRegions]

"""
$(SIGNATURES)
Type of coordinates in grid
"""
coord_type(grid::ExtendableGrid{Tc, Ti}) where {Tc,Ti}=Tc

"""
$(SIGNATURES)
Type of indices
"""
index_type(grid::ExtendableGrid{Tc, Ti}) where {Tc,Ti}=Ti

"""
$(SIGNATURES)

Map a function onto node coordinates of grid
"""
function map(f::Function, grid::ExtendableGrid)
    coord=grid[Coordinates]
    dim=dim_space(grid)
    if dim==1
        @views map(f,coord[1,:])
    elseif dim==2
        @views map(f,coord[1,:], coord[2,:])
    else
        @views map(f,coord[1,:], coord[2,:], coord[3,:])
    end
end



