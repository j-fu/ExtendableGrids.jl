##################################################################
# Abstract types

"""
$(TYPEDEF)
Apex type for grid components.
"""
abstract type AbstractGridComponent <: AbstractExtendableGridApexType end


"""
$(TYPEDEF)
1D Array of floating point data
"""
abstract type AbstractGridFloatArray1D <: AbstractGridComponent end

"""
$(TYPEDEF)
2D Array of floating point data
"""
abstract type AbstractGridFloatArray2D <: AbstractGridComponent end

"""
$(TYPEDEF)
1D Array of interger data
"""
abstract type AbstractGridIntegerArray1D <: AbstractGridComponent end

"""
$(TYPEDEF)
2D Array of integer data
"""
abstract type AbstractGridIntegerArray2D <: AbstractGridComponent end

"""
$(TYPEDEF)
Integer number
"""
abstract type AbstractGridIntegerConstant <: AbstractGridComponent end

"""
$(TYPEDEF)
Floating point number
"""
abstract type AbstractGridFloatConstant <: AbstractGridComponent end


"""
$(TYPEDEF)

Any kind of adjacency between grid components
"""
abstract type AbstractGridAdjacency <: AbstractGridComponent end


"""
$(TYPEDEF)

Array of element geometry information. 
"""
abstract type AbstractElementGeometries <: AbstractGridComponent end


"""
$(TYPEDEF)

Array of element region number information. 
"""
abstract type AbstractElementRegions <: AbstractGridComponent end


##################################################################
# Component key types

"""
$(TYPEDEF)

Node coordinates
"""
abstract type Coordinates <: AbstractGridFloatArray2D end

"""
$(TYPEDEF)

Adjacency describing nodes per grid cell
"""
abstract type CellNodes <: AbstractGridAdjacency end

"""
$(TYPEDEF)

Adjacency describing nodes per grid boundary face
"""
abstract type BFaceNodes <: AbstractGridAdjacency end


"""
$(TYPEDEF)

Description of cell geometries
"""
abstract type CellGeometries <: AbstractElementGeometries end

"""
Description of boundary face geometries

$(TYPEDEF)
"""
abstract type BFaceGeometries <: AbstractElementGeometries end

"""
$(TYPEDEF)

Cell region number per cell
"""
abstract type CellRegions <: AbstractElementRegions end

"""
$(TYPEDEF)

Boundary region number per boundary face
"""
abstract type BFaceRegions <: AbstractElementRegions end

"""
$(TYPEDEF)

Number of cell regions
"""
abstract type NumCellRegions <: AbstractGridIntegerConstant end

"""
$(TYPEDEF)

Number of boundary face regions 
"""
abstract type NumBFaceRegions <: AbstractGridIntegerConstant end

"""
$(TYPEDEF)

Coordinate system
"""
abstract type CoordinateSystem <: AbstractGridComponent end

############################################################
# Grid type

"""
$(TYPEDEF)
Grid type wrapping Dict
"""
mutable struct ExtendableGrid{Tc,Ti}
    components::Dict{Type{<:AbstractGridComponent},Any}
    ExtendableGrid{Tc,Ti}() where{Tc,Ti} =new(Dict{Type{<:AbstractGridComponent},Any}())
end


#############################################################
# 
"""
````
const ElementInfo{T}=Union{Vector{T},VectorOfConstants{T}}`
````

Union type for element information arrays. If all elements have
the same information, it can be stored in an economical form
as a [`VectorOfConstants`](@ref).
"""
const ElementInfo{T}=Union{Vector{T},VectorOfConstants{T}}


############################################################
# Generic component methods

"""
$(TYPEDSIGNATURES)
Default veryform method.

"veryform"  means "verify and/or transform"  and is called to check
and possibly transform components to be added to the grid via `setindex!`.

The default method just passes data through.
"""
veryform(grid::ExtendableGrid,v,::Type{<:AbstractGridComponent})=v


"""
$(TYPEDSIGNATURES)

Set new grid component
"""
Base.setindex!(grid::ExtendableGrid,v,T::Type{<:AbstractGridComponent})= grid.components[T]=veryform(grid,v,T)


"""
$(TYPEDSIGNATURES)

To be called by getindex. This triggers lazy creation of 
non-existing gridcomponents
"""
Base.get!(grid::ExtendableGrid,T::Type{<:AbstractGridComponent})= get!( ()->veryform(grid,instantiate(grid,T),T),  grid.components,T)


"""
$(TYPEDSIGNATURES)

"Hook" for methods instantiating lazy components. 
"""
function instantiate end


"""
````
Base.getindex(grid::ExtendableGrid,T::Type{<:AbstractGridComponent})
````

Generic method for obtaining grid component.

This method is mutating in the sense that non-existing grid components
are created on demand.

Due to the fact that components are stored as Any the return
value triggers type instability.
"""
Base.getindex(grid::ExtendableGrid,T::Type{<:AbstractGridComponent})=get!(grid,T)


############################################################
# Verfication of inserted data

"""
````
veryform(grid::ExtendableGrid{Tc,Ti},v,T::Type{<:AbstractGridAdjacency}) where{Tc,Ti}
````

Check proper type of adjacencies upon insertion
"""
veryform(grid::ExtendableGrid{Tc,Ti},v,T::Type{<:AbstractGridAdjacency}) where{Tc,Ti}= typeof(v)<:Adjacency{Ti} ? v : throw("Type mismatch")


############################################################
# Type stable component access

"""
````
Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridAdjacency})::Adjacency{Ti} where{Tc,Ti}
````

Type stable return of adjacency component
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridAdjacency})::Adjacency{Ti} where{Tc,Ti}
    get!(grid,T)
end

"""
````
Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridFloatArray2D})::Array{Tc,2} 
````

Type stable c method to obtain 2D array from grid
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridFloatArray2D})::Array{Tc,2} where{Tc,Ti}
    get!(grid,T)
end

"""
````
Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridFloatArray1D})::Array{Tc,1} where{Tc,Ti}
````

Type stable method to obtain 1D array from grid
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridFloatArray1D})::Array{Tc,1} where{Tc,Ti}
    get!(grid,T)
end


"""
````
Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractElementGeometries})::ElementInfo{DataType}
````

Type stable method to obtain  element type
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractElementGeometries})::ElementInfo{DataType} where{Tc,Ti}
    get!(grid,T)
end

"""
````
Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractElementRegions})::ElementInfo{Ti} where{Tc,Ti}
````

Type stable method to obtain element region number
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractElementRegions})::ElementInfo{Ti} where{Tc,Ti}
    get!(grid,T)
end


"""
````
Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridIntegerConstant})::Ti where{Tc,Ti}
````

Type stable method to obtain  integer constant stored on elements
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridIntegerConstant})::Ti where{Tc,Ti}
    get!(grid,T)
end

"""
````
Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridFloatConstant})::Tc where{Tc,Ti}
````

Type stable  method to obtain float constant stored on elements
"""
function Base.getindex(grid::ExtendableGrid{Tc,Ti},T::Type{<:AbstractGridFloatConstant})::Tc where{Tc,Ti}
    get!(grid,T)
end


############################################################
# Instantiation methods

"""
$(TYPEDSIGNATURES)

Instantiate number of cell regions
"""
instantiate(grid, ::Type{NumCellRegions})=maximum(grid[CellRegions])

"""
$(TYPEDSIGNATURES)

Instantiate number of bface regions
"""
instantiate(grid, ::Type{NumBFaceRegions})=maximum(grid[BFaceRegions])



#############################################################
# General methods
"""
$(TYPEDSIGNATURES)

Print the hierarchy of grid component key types (subtypes of [`AbstractGridComponent`](@ref). 
This includes additionally user defined subptypes.
"""
gridcomponents()=AbstractTrees.print_tree(AbstractGridComponent,5,indicate_truncation=false)



"""
$(TYPEDSIGNATURES)

Keys in grid
"""
Base.keys(g::ExtendableGrid)=Base.keys(g.components)


"""
$(TYPEDSIGNATURES)

Check if key is in grid
"""
Base.haskey(g::ExtendableGrid,k) = Base.haskey(g.components,k)


"""
$(SIGNATURES)

Space dimension of grid
"""
dim_space(grid::ExtendableGrid)= size(grid[Coordinates],1)


"""
$(SIGNATURES)

Grid dimension dimension of grid (larges element dimension)
"""
dim_grid(grid::ExtendableGrid)=  dim_element(grid[CellGeometries][1])


"""
$(SIGNATURES)

Number of nodes in grid
"""
num_nodes(grid::ExtendableGrid)= size(grid[Coordinates],2)


"""
$(TYPEDSIGNATURES)

Number of cells in grid
"""
num_cells(grid::ExtendableGrid)= num_sources(grid[CellNodes])


"""
$(TYPEDSIGNATURES)

Number of boundary faces in grid.
"""
num_bfaces(grid::ExtendableGrid)= haskey(grid,BFaceNodes) ? num_sources(grid[BFaceNodes]) : 0


"""
$(TYPEDSIGNATURES)

Maximum  cell  region number
"""
num_cellregions(grid::ExtendableGrid)=grid[NumCellRegions]


"""
$(TYPEDSIGNATURES)

Maximum  boundary face region numbers
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
$(TYPEDSIGNATURES)

Map a function onto node coordinates of grid
"""
function Base.map(f::Function, grid::ExtendableGrid)
    coord=grid[Coordinates]
    dim=dim_space(grid)
    if dim==1
        @views Base.map(f,coord[1,:])
    elseif dim==2
        @views Base.map(f,coord[1,:], coord[2,:])
    else
        @views Base.map(f,coord[1,:], coord[2,:], coord[3,:])
    end
end


function Base.show(io::IO,grid::ExtendableGrid)
    if num_edges(grid)>0
        str=@sprintf("%s;\ndim: %d nodes: %d cells: %d bfaces: %d, edges: %d\n",
                     typeof(grid),dim_space(grid),num_nodes(grid), num_cells(grid), num_bfaces(grid), num_edges(grid))
    else
        str=@sprintf("%s;\ndim: %d nodes: %d cells: %d bfaces: %d\n",
                     typeof(grid),dim_space(grid),num_nodes(grid), num_cells(grid), num_bfaces(grid))
    end
    println(io,str)
end

"""
$(SIGNATURES)

Recursively check seeming equality of two grids. Seemingly means 
that long arrays are only compared via random samples
"""
function seemingly_equal(grid1::ExtendableGrid, grid2::ExtendableGrid)
    for key in keys(grid1)
        if !haskey(grid2,key)
            return false
        end
        if !seemingly_equal(grid1[key],grid2[key])
            return false
        end
    end
    return true
end

"""
$(SIGNATURES)

Check for seeming equality of two arrays by random sample.
"""
function seemingly_equal(array1::AbstractArray, array2::AbstractArray)
    if size(array1)!=size(array2)
        return false
    end
    l=length(array1)
    ntests=Float64(min(l,50))
    p=min(0.5,ntests/l)
    for i in randsubseq( 1:l,p)
        if !seemingly_equal(array1[i],array2[i])
            return false
        end
    end
    true
end

seemingly_equal(a1::VariableTargetAdjacency, a2::VariableTargetAdjacency)=seemingly_equal(a1.colentries,a2.colentries)&& seemingly_equal(a1.colstaert, a2.colstart)
seemingly_equal(x1::Type, x2::Type)=(x1==x2)
seemingly_equal(x1::Number, x2::Number)=(x1≈x2)
seemingly_equal(x1::Any, x2::Any)=(x1==x2)


Base.extrema(grid::ExtendableGrid)=Base.extrema(grid[Coordinates])

function bbox(grid)
    e=extrema(grid)
    map(a->a[1],e),map(a->a[2],e)
end
