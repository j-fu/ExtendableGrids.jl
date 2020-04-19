"""
# Extendable  containers with type stable content access and lazy content creation via the Julia type system.

## Problem to be addressed

In certain contexts it is desirable to use containers with core components
which are user extendable and allow for type stable component acces. Moreover,
some components are necessary on demand only, so they should be created lazily.
Furthermore, there should be a kind of safety protocol which prevents errors
from typos in component names etc.

Julia default data structures do not provide these properties.

### `struct` 
  - Julia structs with proper field type annotations guarantee type stability
  - Julia structs are not extendable, fields and their types are fixed upon definition
  - If we don't fix types of struct fields they become Any and a source 
    for type instability
  - The situation could be fixed if `getfield` could be overloaded but it cant't

### `Dict`
  - Plain Dicts with flexible value types are a source of type instability
  - Typical use of Dicts with strings as keys needs a meta protocol to handle
    semantics of keys which at the end probably hinges on string comparison which
    will make things slow
  - Same for the implementation of a lazy evaluation protocol

## Proposed solution:

Harness the power of the Julia type system: 
- Use a struct containing a  Dict with DataType as keys. Every key is a type.
- Use the type system to dispatch between  `getindex`/`setindex!` methods for keys
- Use type hierarchies to manage different different value classes
- Extension requires declaring new types, keys can be only existing types almost removing
  typos as sources for errors
- Lazy extension is managed bye an  `instantiate` method called by `getindex` if necessary
- Component access is made type stable by type dependent `getindex` methods
- Component insertion is made safe by having  `setindex!`  calling a `veryform` method

### Pros
See above ...

### Cons
- Implemented using a Dict, so access is inherently slower than access to a component
  of a struct. Therefore it is not well suited for inner loops.
    
"""


"""
Apex type for grid components
"""
abstract type AbstractGridComponent end

"""
2D Array on grid components
"""
abstract type AbstractGridArray2D <: AbstractGridComponent end

"""
1D Array on grid components
"""
abstract type AbstractGridArray1D <: AbstractGridComponent end

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
Basic Grid components with classification
"""
abstract type Coordinates <: AbstractGridArray2D end

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
    components::Dict{DataType,Any}
    ExtendableGrid{Tc,Ti}() where{Tc,Ti} =new(Dict{AbstractGridComponent,Any}())
end

"""
Default veryform method.

"veryform"  means "verify and/or transform"  and is called to check
and possibly transform components to be added to the grid
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
"Hook" for method instantiating lazy components. 
We need to have at least one method here to define a function which then can 
get more user-defined methodsa
"""
instantiate(grid, ::Type{T}) where T = nothing

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
Type specific method to obtain grid vector
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{AbstractGridArray2D})::Array{Tc,2} where{Tc,Ti}
    get!(grid,T)
end

"""
Type specific method to obtain grid vector
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{AbstractGridArray1D})::Array{Tc,1} where{Tc,Ti}
    get!(grid,T)
end

"""
Type specific method to obtain  element type
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{AbstractElementTypes})::ElementInfo{AbstractElementType} where{Tc,Ti}
    get!(grid,T)
end

"""
Type specific method to obtain  element type
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{AbstractElementRegions})::ElementInfo{Ti} where{Tc,Ti}
    get!(grid,T)
end



