"""
$(TYPEDEF)

Grid component key type for storing node parents (=ids of nodes in ParentGrid) in an array
"""
abstract type NodeParents <: AbstractGridIntegerArray1D end

"""
$(TYPEDEF)

Grid component key type for storing parent faces
(only for SubGrid relation when FaceNodes is instantiated)
"""
abstract type FaceParents <: AbstractGridIntegerArray1D end

"""
$(TYPEDEF)

Grid component key type for storing parent bfaces
"""
abstract type BFaceParents <: AbstractGridIntegerArray1D end

"""
$(TYPEDEF)

Grid component key type for storing parent cells
"""
abstract type CellParents <: AbstractGridIntegerArray1D end

"""
$(TYPEDEF)

Grid component key type for storing parent grid
"""
abstract type ParentGrid <: AbstractGridComponent end

"""
$(TYPEDEF)

Grid component key type for storing parent grid relationship
"""
abstract type ParentGridRelation <: AbstractGridComponent end

"""
$(TYPEDEF)

Grid component key type for indicating that grid is a subgrid of the parentgrid
"""
abstract type SubGrid <: ParentGridRelation end

"""
$(TYPEDEF)

Grid component key type for indicating that grid is a boundary subgrid of the parentgrid
"""
abstract type BoundarySubGrid <: ParentGridRelation end

"""
$(TYPEDEF)

Grid component key type for indicating that grid is a refinement of the parentgrid
"""
abstract type RefinedGrid <: ParentGridRelation end


"""
$(TYPEDSIGNATURES)

Default transform for subgrid creation
"""
function _copytransform!(a::AbstractArray,b::AbstractArray)
    for i=1:length(a)
        a[i]=b[i]
    end
end

struct XIPair{Tv, Ti}
    x::Tv
    i::Ti
end

# Comparison method for sorting
Base.isless(x::XIPair, y::XIPair) = (x.x < y.x)



"""
$(TYPEDSIGNATURES)

Create subgrid from list of regions.

- `parent`: parent grid 
- `subregions`:  Array of subregions
- `transform` (kw parameter): transformation function between
   grid and subgrid coordinates acting on one point.
   Default: `copytransform`
- `boundary`: if true, create codimension 1 subgrid from boundary region.
- `project`: project coordinates onto  subgrid dimension

A subgrid is of type `ExtendableGrid` and stores two additional components:
[`ParentGrid`](@ref) and [`NodeParents`](@ref)

"""
function subgrid(parent,
                 subregions::AbstractArray;
                 transform::Function=_copytransform!,
                 boundary=false,
                 project=true)

    Tc=coord_type(parent)
    Ti=index_type(parent)

    #
    # TODO: make a flag array here
    #
    @inline function insubregions(xreg)
        for i in eachindex(subregions)
            if subregions[i]==xreg
                return true
            end
        end
        return false
    end

    
    if boundary
        xregions=parent[BFaceRegions]
        xnodes=parent[BFaceNodes]
        sub_gdim=dim_grid(parent)-1
        xct=parent[BFaceGeometries]
        sub_gdim=dim_grid(parent)-1
    else
        xregions=parent[CellRegions]
        xnodes=parent[CellNodes]
        xct=parent[CellGeometries]
        sub_gdim=dim_grid(parent)
    end
    
    nodemark=zeros(Ti,num_nodes(parent))
    
    nsubcells=0
    nsubnodes=0
    cellparents=zeros(Ti,0)
    for icell in eachindex(xregions)
        if insubregions(xregions[icell])
            nsubcells+=1
            for inode=1:num_targets(xnodes,icell)
                ipnode=xnodes[inode,icell]
                if nodemark[ipnode]==0
                    nsubnodes+=1
                    nodemark[ipnode]=nsubnodes
                end
            end
            push!(cellparents, icell)
        end
    end
    
    sub_xnodes=VariableTargetAdjacency(Ti)
    sub_nip=zeros(Ti,nsubnodes)
    sub_ct=Vector{ElementGeometries}(undef,0)
    sub_cr=Vector{Ti}(undef,0)
    for inode in eachindex(nodemark)
        if nodemark[inode]>0
            sub_nip[nodemark[inode]]=inode
        end
    end
    
    isubcell=0
    for icell in eachindex(xregions)
        if insubregions(xregions[icell])
            ncn=num_targets(xnodes,icell)
            col=zeros(Ti,0)
            for inode=1:ncn
                push!(col,nodemark[xnodes[inode,icell]])
            end
            append!(sub_xnodes,col)
            push!(sub_ct,xct[icell])
            push!(sub_cr,xregions[icell])
        end
    end

    if project
        sub_coord=zeros(Tc,sub_gdim,nsubnodes)
    else
        sub_coord=zeros(Tc,dim_space(parent),nsubnodes)
    end

    coord=parent[Coordinates]
    @views for inode=1:nsubnodes
        transform(sub_coord[:,inode],coord[:,sub_nip[inode]])
    end
    subgrid=ExtendableGrid{Tc,Ti}()
    subgrid[Coordinates]=sub_coord
    subgrid[CellRegions]=sub_cr
    subgrid[CellGeometries]=sub_ct
    subgrid[CellNodes]=tryfix(sub_xnodes)
    subgrid[ParentGrid]=parent
    subgrid[NodeParents]=sub_nip
    subgrid[CellParents]=cellparents
    subgrid[ParentGridRelation]=boundary ? BoundarySubGrid : SubGrid

    if boundary
        subgrid[NumBFaceRegions]=0
        subgrid[BFaceRegions]=Ti[]
        subgrid[BFaceGeometries]=ElementGeometries[]
        subgrid[BFaceNodes]=Matrix{Ti}(undef,sub_gdim,0)
        subgrid[NumBFaceRegions]=0
    else
        bfacenodes=parent[BFaceNodes]
        bfaceregions=parent[BFaceRegions]
        bfacetypes=parent[BFaceGeometries]
        bfacecells=parent[BFaceCells]
        bfaceparents=zeros(Ti,0)
        
        sub_bfacenodes=VariableTargetAdjacency(Ti)
        sub_bfaceregions=Vector{Ti}(undef,0)
        sub_bfacetypes=Vector{ElementGeometries}(undef,0)
        
        for ibface in eachindex(bfaceregions)
            nbn=num_targets(bfacenodes,ibface)
            insubgrid=true
            for inode=1:nbn
                if nodemark[bfacenodes[inode,ibface]]==0
                    insubgrid=false
                    continue
                end
            end
            # All nodes may be in subgrid, but this still doesn't mean
            # yet that the bface is in subgrid - we need to check if the
            # neigboring cells are in the subgrid
            # TODO: this even may be sufficient!
            insubgrid=false
            for itarget=1:num_targets(bfacecells,ibface)
                icell=bfacecells[itarget, ibface]
                if insubregions(xregions[icell])
                    insubgrid=true
                    continue
                end
            end
            if insubgrid
                col=zeros(Ti,0)
                for inode=1:nbn
                    push!(col,nodemark[bfacenodes[inode,ibface]])
                end
                append!(sub_bfacenodes,col)
                push!(sub_bfacetypes,bfacetypes[ibface])
                push!(sub_bfaceregions,bfaceregions[ibface])
                push!(bfaceparents, ibface)
            end
        end
        subgrid[BFaceParents] = bfaceparents
        subgrid[BFaceRegions]=sub_bfaceregions
        subgrid[BFaceGeometries]=sub_bfacetypes
        if length(sub_bfaceregions) > 1
            subgrid[BFaceNodes]=tryfix(sub_bfacenodes)
            subgrid[NumBFaceRegions]=maximum(sub_bfaceregions)
        else
            subgrid[BFaceNodes]=zeros(Ti, 2, 0)
            subgrid[NumBFaceRegions]=0
        end
    end
    subgrid[CoordinateSystem]=parent[CoordinateSystem]

    if sub_gdim == 1
        # Sort nodes of grid for easy plotting
        X=view(subgrid[Coordinates],1,:)
        nx=length(X)
        I=subgrid[NodeParents]
        xipairs=[XIPair{Tc,Ti}(X[i],I[i]) for i=1:nx]
        sort!(xipairs, 1,nx, Base.QuickSort, Base.Forward)
        for i=1:nx
            X[i]=xipairs[i].x
            I[i]=xipairs[i].i
        end
    end
    
    subgrid
end




"""
$(TYPEDEF)
Vector view on subgrid

$(TYPEDFIELDS)
"""
struct SubgridVectorView{Tv,Ti} <: AbstractVector{Tv}
    sysarray::AbstractVector{Tv}
    node_in_parent::Vector{Ti}
end

##################################################################
"""
$(TYPEDSIGNATURES)

Create a view of the vector on a subgrid.
"""
Base.view(a::AbstractVector,subgrid::ExtendableGrid)  = SubgridVectorView(a,subgrid[NodeParents])


##############################################################################
"""
$(TYPEDSIGNATURES)

Accessor method for subgrid vector view.
"""
Base.getindex(aview::SubgridVectorView,inode::Integer) = aview.sysarray[aview.node_in_parent[inode]]

##############################################################################
"""
$(TYPEDSIGNATURES)

Accessor method for subgrid vector view.
"""
@inline function Base.setindex!(aview::SubgridVectorView,v,inode::Integer)
    aview.sysarray[aview.node_in_parent[inode]]=v
    return aview
end

##################################################################
"""
$(TYPEDSIGNATURES)
    
Return size of vector view.
"""
Base.size(a::SubgridVectorView)=(size(a.node_in_parent,1),)



