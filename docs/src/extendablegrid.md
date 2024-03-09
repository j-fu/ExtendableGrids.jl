# Extendable grid


An ExtendableGrid in form of a dictionary with types as keys and type stable value access.
This means that grid components are accessed as dict entries, e.g. `grid[Coordinates]`.
The rationale behind this decision is described [below](@ref TDict).

```@contents
Pages = ["extendablegrid.md"]
Depth = 2:6
```

## Extendable grid notations
A grid is assumed to be a subset of components of a polyhedral complex in d-dimensional space.
We distinguish the following element classes characterized by their dimension:

| Element class | Meaning                                                                    |
| :--           | :--                                                                        |
| Node          | 0-dimensional node                                                         |
| Edge          | 1-dimensional line connecting two neigboring nodes                         |
| Face          | codimension 1 object separating a cell from outer space or neigboring cell |
| Cell          | codimension 0 object                                                       |
| BFace         | Face situated at inner or domain boundary                                  |
| Region        | number to be used to characterize subdomains, contacts etc.                |
 
## Grid components

Grid components are accessed like Dict entries, the keys must be subtypes of [`AbstractGridComponent`](@ref).

### Basic set of grid components
Upon construction, an ExtendableGrid needs to be provided with the basic set of grid components denoted by
the following component type keys:

| Component type key          | Meaning                                                                                              |
| :---------------            | :--------------------------------------------------------------------------------------------------- |
| [`Coordinates`](@ref)       | Coordinates of the vertices of the grid cells                                                        |
| [`CellNodes`](@ref)         | Adjacency describing the nodes of grid cell                                                          |
| [`CellGeometries`](@ref)    | Abstract array of subtypes of AbstractElementGeometry describing the geometry of each cell           |
| [`CellRegions`](@ref)       | Abstract array of integers describing region numbers                                                 |
| [`BFaceNodes`](@ref)        | Adjacency structure describing the nodes corresponding to each grid boundary face                    |
| [`BFaceGeometries`](@ref)   | Abstract array of subtypes of AbstractElementGeometry describing the geometry of each boundary face  |
| [`BFaceRegions`](@ref)      | Abstract array of integers describig region numbers                                                  |
| [`CoordinateSystem `](@ref) | Abstract type describing the coordinate system to be used                                            |



### Hierarchy of component type keys

The list of components can be printed using the [`gridcomponents`](@ref) method.
```@example
using ExtendableGrids # hide
gridcomponents() #hide
```

### Additional components
Additional components can be added by defining  a subtype of [`AbstractGridComponent`](@ref) or
a fitting subtype thereof, and assigning the value to the corresponding Dict entry:

```@example
using ExtendableGrids # hide
g=simplexgrid([1,2,3,4.0])
abstract type MyComponent <: AbstractGridComponent end
g[MyComponent]=13
show(g)
```

Alternatively, component creation can be perfomed lazily. For this
purpose one needs to define an `instantiate` method:


```@example
using ExtendableGrids # hide
abstract type NodeCells <: AbstractGridAdjacency end
ExtendableGrids.instantiate(grid, ::Type{NodeCells})=atranspose(grid[CellNodes])
g=simplexgrid([1,2,3,4.0])
show(g[NodeCells])
```


## Grid API

```@autodocs
Modules = [ExtendableGrids]
Pages = ["extendablegrid.jl"]
```



## [The TDict interface pattern](@id TDict)

Here we describe the idea behind the data structure used in this package.
TDict means: extendable containers with type stable content access and lazy content creation via the Julia type system.

### Problem to be addressed

In certain contexts it is desirable to use containers with core components
which are user extendable and allow for type stable component acces. Moreover,
some components are necessary on demand only, so they should be created lazily.
Furthermore, there should be a kind of safety protocol which prevents errors
from typos in component names etc.

Julia default data structures do not provide these properties.

#### `struct` 
  - Julia structs with proper field type annotations guarantee type stability
  - Julia structs are not extendable, fields and their types are fixed upon definition
  - If we don't fix types of struct fields they become Any and a source 
    for type instability
  - The situation could be fixed if `getfield` could be overloaded but it cant't

#### `Dict`
  - Plain Dicts with flexible value types are a source of type instability
  - Dicts with strings as keys needs a meta protocol to handle
    semantics of keys which at the end probably hinges on string comparison which
    will make things slow
  - Dicts with symbols as keys still need this meta protocol
  - Same for the implementation of a lazy evaluation protocol
  - If a dict contains components of different types, component access will not be typestable

### Proposed solution:

Harness the power of the Julia type system: 
- Use a struct containing a  Dict with DataType as keys. Every key is a type.
- Use type hierarchies to manage different  value classes
- Use the type system to dispatch between  `getindex`/`setindex!` methods for keys
- Extension requires declaring new types, keys can be only existing types almost removing
  typos as sources for errors
- Lazy extension is managed bye an  `instantiate` method called by `getindex` if necessary
- Component access is made type stable by type dispatched`getindex` methods
- Component insertion is made safe by having  `setindex!`  calling a `veryform` method

#### Pros
See above ...

#### Cons
- Implemented using a Dict, so access is inherently slower than access to a component
  of a struct. Therefore it is not well suited for inner loops.

