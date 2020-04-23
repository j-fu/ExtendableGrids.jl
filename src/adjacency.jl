
##########################################################################################################
"""
   Adjacency interface.

   This handles adjacency matrices between entities of polyhedral complexes, e.g.
   nodes, cells, edges etc.
     
   An adjacency is described by an Adjacency matrix, which is a sparse matrix
   whose entries a 0 or 1. While such a matrix always can be stored
   as a SparseMatrixCSC, in general this would be a waste of storage.
   
   For the general case, it is sufficient to only store the column start
   indieces and the column entries (row numbers), and to implicitely assume
   that nonzero entries are 1. This kind of storage is realised in a
   VariableTargetAdjacency.

   In many cases, this can be compressed even more, if each column has the
   same length. In that case, a Matrix is sufficient to store the data.
   This is the usual base for implementing FEM/FVM assembly, and the interface
   for the general case should be similar.

   From these ideas we develop the following interface for an adjacency a.

   In order to avoid name confusion, we introduce the following notation which 
   should be consistent with the use in assembly loops.

   source:  source of adjacency link
   target:  target of adjacency link
 
   E.g. the cell-node adjacency for FEM assembly links  a number of
   cells with a collection of nodes.  The cells are the sources,
   and the targets are the nodes. 

   getindex(a,i,isource) aka a[i,isource]: return i-th target of  source j
   num_sources(a): overall number of sources, e.g. number of cells
   num_targets(a): overall number of targets
   num_targets(a,isource): number of targets for source given by isource
   num_links(a): number of links aka nonzero entries of adjacency matrix
   show(a): print stuff

   Further API ideas:
   - Convert between Matrix and Variable target stuff using 0 entries as "padding"
"""
struct VariableTargetAdjacency{T}
    colentries::Vector{T}
    colstart::Vector{T}
end

function Base.:(==)(a::VariableTargetAdjacency{Ta}, b::VariableTargetAdjacency{Tb}) where {Ta,Tb}
    Ta==Tb && a.colentries==b.colentries &&  a.colstart==b.colstart
end             

"""
Create an empty VariableTargetAdjacency
"""
VariableTargetAdjacency(t::Type{T}) where T=VariableTargetAdjacency{T}(Vector{T}(undef,0),[one(T)])

"""
Create an empty VariableTargetAdjacency with default type
"""
VariableTargetAdjacency()=VariableTargetAdjacency(Int64)


"""
    Create a VariableTargetAdjacency from Matrix
"""
VariableTargetAdjacency(m::Matrix{T}) where T=VariableTargetAdjacency{T}(vec(m),collect(1:size(m,1):size(m,1)*size(m,2)+1))

"""
    Show adjacency (in trasposed form; preliminary)
"""
function Base.show(io::IO,adj::VariableTargetAdjacency)
    for isource=1:num_sources(adj)
        for itarget=1:num_targets(adj,isource)
            print(adj[itarget,isource], " ")
        end
        println()
    end
end

"""
    Access adjacency as if it is a 2D Array
"""
Base.getindex(adj::VariableTargetAdjacency,i,isource)=adj.colentries[adj.colstart[isource]+i-1]

Base.getindex(adj::VariableTargetAdjacency,::Colon,isource)=adj.colentries[adj.colstart[isource]:adj.colstart[isource+1]-1]

"""
    Number of targets for given source
"""
num_targets(adj::VariableTargetAdjacency,isource)=adj.colstart[isource+1]-adj.colstart[isource]

"""
    Number of sources in adjacency
"""
num_sources(adj::VariableTargetAdjacency)=length(adj.colstart)-1

"""
    Number of targeta
"""
num_targets(adj::VariableTargetAdjacency)=maximum(adj.colentries)

"""
    Number of links
"""
num_links(adj::VariableTargetAdjacency)=length(adj.colentries)

"""
    Maximum number of targets per source
"""
max_num_targets_per_source(adj::VariableTargetAdjacency)=maximum(adj.colstart[2:end].-adj.colstart[1:end-1])


"""
    Append a column to adjacency.
"""
function Base.append!(adj::VariableTargetAdjacency,column)
    for i=1:length(column)
        push!(adj.colentries,column[i])
    end
    push!(adj.colstart,length(adj.colentries)+1)
end

"""
    Use Matrix to store fixed target adjacency
"""
const FixedTargetAdjacency=Matrix

"""
    Number of targets per source if adjacency is a matrix
"""
num_targets(adj::FixedTargetAdjacency,isource)=size(adj)[1]

"""
    Number of sources in adjacency
"""
num_sources(adj::FixedTargetAdjacency)=size(adj)[2]

"""
    Overall number of targets 
"""
num_targets(adj::FixedTargetAdjacency)=maximum(vec(adj))

"""
    Number of entries
"""
num_links(adj::FixedTargetAdjacency)=length(adj)

"""
    Maximum number of targets per source
"""
max_num_targets_per_source(adj::FixedTargetAdjacency)=size(adj,1)


const Adjacency{T}=Union{FixedTargetAdjacency{T},VariableTargetAdjacency{T}}

Adjacency{T}(a::FixedTargetAdjacency{T}) where T =a
Adjacency{T}(a::VariableTargetAdjacency{T}) where T =a

function atranspose(adj::Adjacency{T}) where T
    # 0th pass: calculate number of rows !!! todo: how to call ?
    t_adj=VariableTargetAdjacency(zeros(T,num_links(adj)),zeros(T,num_targets(adj)+1))
    
    # 1st pass: calculate new column sizes, store them in t_adj.colstart
    for isource=1:num_sources(adj)
        for itarget=1:num_targets(adj,isource)
            t_adj.colstart[adj[itarget,isource]]+=1
        end
    end
    
    # 2nd pass: calculate initial assembly
    # addresses in adj.ja by summing up deltas
    # store them a the ends of  the columns:
    # in t_adj.ja[t_adj.colstart[isource+1]-1].
    # it will be overwritten when the last element of
    # the column is assembled.
    # We get automatically increasing  column indices.
    
    delta=t_adj.colstart[1];
    t_adj.colstart[1]=1;
    for isource=2:num_sources(t_adj)
	save=t_adj.colstart[isource];
	t_adj.colstart[isource]=t_adj.colstart[isource-1]+delta; 
	if delta>0
            t_adj.colentries[t_adj.colstart[isource]-1]=t_adj.colstart[isource-1];
        end
	delta=save;
    end
    isource=num_sources(t_adj)+1
    t_adj.colstart[isource]=t_adj.colstart[isource-1]+delta;
    if delta>0
        t_adj.colentries[t_adj.colstart[isource]-1]=t_adj.colstart[isource-1];
    end
    
    # 3rd pass: assemble new columns
    for isource=1:num_sources(adj)
        for itarget=1:num_targets(adj,isource)
	    asm_idx=t_adj.colstart[adj[itarget,isource]+1]-1;
	    asm_loc=t_adj.colentries[asm_idx];
	    t_adj.colentries[asm_idx]+=1; 
	    t_adj.colentries[asm_loc]=isource;
        end
    end
    return t_adj
end

