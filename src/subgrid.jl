"""
$(TYPEDEF)

Grid component key type for storing node in parent array
"""
abstract type NodeInParent <: AbstractGridIntegerArray1D end

"""
$(TYPEDEF)

Grid component key type for storing parent grid
"""
abstract type ParentGrid <: AbstractGridComponent end

"""
$(TYPEDSIGNATURES)

Default transform for subgrid creation
"""
function _copytransform!(a::AbstractArray,b::AbstractArray)
    for i=1:length(a)
        a[i]=b[i]
    end
end

"""
$(TYPEDSIGNATURES)

Create subgrid from list of regions.

- `parent`: parent grid 
- `subregions`:  Array of subregions
- `transform` (kw parameter): transformation function between
   grid and subgrid coordinates acting on one point.
   Default: `copytransform`
- `boundary`: if true, create codimension 1 subgrid from boundary region.

A subgrid is of type `ExtendableGrid` and stores two additional components:
[`ParentGrid`](@ref) and [`NodeInParent`](@ref)

"""
function subgrid(parent,
                 subregions::AbstractArray;
                 transform::Function=_copytransform!,
                 boundary=false,
                 project=true)

    Tc=coord_type(parent)
    Ti=index_type(parent)
    
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
    subgrid[NodeInParent]=sub_nip

    if boundary
        subgrid[NumBFaceRegions]=0
    else
        bfacenodes=parent[BFaceNodes]
        bfaceregions=parent[BFaceRegions]
        bfacetypes=parent[BFaceGeometries]
        
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
            if insubgrid
                col=zeros(Ti,0)
                for inode=1:nbn
                    push!(col,nodemark[bfacenodes[inode,ibface]])
                end
                append!(sub_bfacenodes,col)
                push!(sub_bfacetypes,bfacetypes[ibface])
                push!(sub_bfaceregions,bfaceregions[ibface])
            end
        end
    
        subgrid[BFaceRegions]=sub_bfaceregions
        subgrid[BFaceGeometries]=sub_bfacetypes
        subgrid[BFaceNodes]=tryfix(sub_bfacenodes)
        subgrid[NumBFaceRegions]=maximum(sub_bfaceregions)
    end
    subgrid[CoordinateSystem]=parent[CoordinateSystem]
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
Base.view(a::AbstractVector,subgrid::ExtendableGrid)  = SubgridVectorView(a,subgrid[NodeInParent])


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



