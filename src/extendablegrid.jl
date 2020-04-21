"""
Apex type for grid components
"""
abstract type AbstractGridComponent end

"""
2D Array on grid components (e.g. coordinates)
"""
abstract type AbstractGridFloatArray2D <: AbstractGridComponent end

"""
1D Array on grid components
"""
abstract type AbstractGridFloatArray1D <: AbstractGridComponent end

"""
Integer number
"""
abstract type AbstractGridIntegerConstant <: AbstractGridComponent end

"""
Floating point  number
"""
abstract type AbstractGridFloatConstant <: AbstractGridComponent end


"""
Any kind of adjacency between grid components
"""
abstract type AbstractGridAdjacency <: AbstractGridComponent end


"""
    ElementInfo{AbstractCellType} on arrays
"""
abstract type AbstractElementTypes <: AbstractGridComponent end

"""
    ElementInfo{Ti} on arrays
"""
abstract type AbstractElementRegions <: AbstractGridComponent end

"""
Basic Grid components with classification via intermediate types
"""
abstract type Coordinates <: AbstractGridFloatArray2D end

abstract type NumCellRegions <: AbstractGridIntegerConstant end
abstract type NumBFaceRegions <: AbstractGridIntegerConstant end

abstract type CellNodes <: AbstractGridAdjacency end
abstract type BFaceNodes <: AbstractGridAdjacency end

abstract type CellTypes <: AbstractElementTypes end
abstract type BFaceTypes <: AbstractElementTypes end

abstract type CellRegions <: AbstractElementRegions end
abstract type BFaceRegions <: AbstractElementRegions end


"""
Grid type wrapping Dict
"""
mutable struct ExtendableGrid{Tc,Ti}
    components::Dict{Type{<:AbstractGridComponent},Any}
    ExtendableGrid{Tc,Ti}() where{Tc,Ti} =new(Dict{Type{<:AbstractGridComponent},Any}())
end

"""
Default veryform method.

"veryform"  means "verify and/or transform"  and is called to check
and possibly transform components to be added to the grid.


The default method just passes data through.
"""
veryform(grid::ExtendableGrid,v,::Type{<:AbstractGridComponent})=v

"""
Check proper type of adjacencies
"""
veryform(grid::ExtendableGrid{Tc,Ti},v,T::Type{<:AbstractGridAdjacency}) where{Tc,Ti}= typeof(v)<:Adjacency{Ti} ? v : throw("Type mismatch")


"""
Set new grid component
"""
Base.setindex!(grid::ExtendableGrid,v,T::Type{<:AbstractGridComponent})= grid.components[T]=veryform(grid,v,T)


"""
To be called by getindex. This triggers lazy creation of 
non-existing gridcomponents
"""
Base.get!(grid::ExtendableGrid,T::Type{<:AbstractGridComponent})= get!( ()->veryform(grid,instantiate(grid,T),T),  grid.components,T)


"""
"Hook" for methods instantiating lazy components. 
See https://white.ucc.asn.au/2020/04/19/Julia-Antipatterns.html
"""
function instantiate end


"""
Generic method for obtaining grid component.

This method is mutating in the sense that non-existing grid components
are created on demand.

Due to the fact that components are stored as Any the return
value triggers type instability.
"""
Base.getindex(grid::ExtendableGrid,T::Type{<:AbstractGridComponent})=get!(grid,T)


"""
Type specific method to obtain adjacency component
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridAdjacency})::Adjacency{Ti} where{Tc,Ti}
    get!(grid,T)
end

"""
Type specific method to obtain 2D array from grid
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridFloatArray2D})::Array{Tc,2} where{Tc,Ti}
    get!(grid,T)
end

"""
Type specific method to obtain 1D array from grid
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridFloatArray1D})::Array{Tc,1} where{Tc,Ti}
    get!(grid,T)
end

"""
Type specific method to obtain  element type
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractElementTypes})::ElementInfo{AbstractElementType} where{Tc,Ti}
    get!(grid,T)
end

"""
Type specific method to obtain  element type
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractElementRegions})::ElementInfo{Ti} where{Tc,Ti}
    get!(grid,T)
end


function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridIntegerConstant})::Ti where{Tc,Ti}
    get!(grid,T)
end

function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridFloatConstant})::Tc where{Tc,Ti}
    get!(grid,T)
end


instantiate(grid, ::Type{NumCellRegions})=maximum(grid[CellRegions])
instantiate(grid, ::Type{NumBFaceRegions})=maximum(grid[BFaceRegions])

######################################################################################################################
#Definition of interface methods for grid.

#  ... die könnten alle als lazyconstants gehen



##########################################################
"""
$(TYPEDSIGNATURES)

Space dimension of grid
"""
dim_space(grid::ExtendableGrid)= size(grid[Coordinates],1)


##########################################################
"""
$(TYPEDSIGNATURES)


Number of nodes in grid
"""
num_nodes(grid::ExtendableGrid)= size(grid[Coordinates],2)


##########################################################
"""
$(TYPEDSIGNATURES)

Number of cells in grid
"""
num_cells(grid::ExtendableGrid)= num_sources(grid[CellNodes])

##########################################################
"""
$(TYPEDSIGNATURES)

Number of edges in grid
"""
num_edges(grid::ExtendableGrid)= num_sources(grid[EdgeNodes])


################################################
"""
$(TYPEDSIGNATURES)

Number of boundary faces in grid.
"""
num_bfaces(grid::ExtendableGrid)= num_sources(grid[BFaceNodes])

################################################
"""
$(TYPEDSIGNATURES)

Maximum  cell  region number
"""
num_cellregions(grid::ExtendableGrid)=grid[NumCellRegions]


################################################
"""
$(TYPEDSIGNATURES)

Maximum  boundary face region number
"""
num_bfaceregions(grid::ExtendableGrid)=grid[NumBFaceRegions]




################################################
"""
$(TYPEDSIGNATURES)

Map a function onto grid coordinates.
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
