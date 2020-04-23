abstract type NodeInParent <: AbstractGridIntegerArray1D end
abstract type ParentGrid <: AbstractGridComponent end

##################################################################
# Default transform for subgrid creation
function _copytransform!(a::AbstractArray,b::AbstractArray)
    for i=1:length(a)
        a[i]=b[i]
    end
end

##################################################################
"""
$(TYPEDSIGNATURES)

Create subgrid of list of regions.
"""
function subgrid(parent,
                 subregions::AbstractArray;
                 transform::Function=_copytransform!,
                 boundary=false)
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
        xct=parent[BFaceTypes]
        sub_gdim=dim_grid(parent)-1
    else
        xregions=parent[CellRegions]
        xnodes=parent[CellNodes]
        xct=parent[CellTypes]
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
    sub_ct=Vector{DataType}(undef,0)
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

    sub_coord=zeros(Tc,sub_gdim,nsubnodes)
    coord=parent[Coordinates]
    @views for inode=1:nsubnodes
        transform(sub_coord[:,inode],coord[:,sub_nip[inode]])
    end
    
    subgrid=ExtendableGrid{Tc,Ti}()
    subgrid[Coordinates]=sub_coord
    subgrid[CellRegions]=sub_cr
    subgrid[CellTypes]=sub_ct
    subgrid[CellNodes]=sub_xnodes
    subgrid[ParentGrid]=parent
    subgrid[NodeInParent]=sub_nip
    subgrid
end




"""
$(TYPEDEF)

Struct holding information for solution array view on subgrid

$(TYPEDFIELDS)
"""
struct XSubgridVectorView{Tv,Ti} <: AbstractVector{Tv}

    sysarray::AbstractVector{Tv}

    node_in_parent::Vector{Ti}
end

##################################################################
"""
$(TYPEDSIGNATURES)

Create a view of the solution array on a subgrid.
"""
Base.view(a::AbstractVector,subgrid::ExtendableGrid)  = XSubgridVectorView(a,subgrid[NodeInParent])


##############################################################################
"""
$(TYPEDSIGNATURES)

Accessor method for subgrid array view.
"""
Base.getindex(aview::XSubgridVectorView,inode::Integer) = aview.sysarray[aview.node_in_parent[inode]]

##############################################################################
"""
$(TYPEDSIGNATURES)

Accessor method for subgrid array view.
"""
@inline function Base.setindex!(aview::XSubgridVectorView,v,inode::Integer)
    aview.sysarray[aview.node_in_parent[inode]]=v
    return aview
end

##################################################################
"""
$(TYPEDSIGNATURES)
    
Return size of solution array view.
"""
Base.size(a::XSubgridVectorView)=(size(a.node_in_parent,1),)



