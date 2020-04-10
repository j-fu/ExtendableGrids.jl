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
   that nonzero entries are 1.

   Im many cases, this can be compressed even more, if each column has the
   same length. In that case, a Matrix is sufficient to store the data.
   This is the usual base for implementing FEM/FVM assembly, and the interface
   for the general case should be similar.

   From these ideas we develop the following interface for an adjacency a.

   In order to avoid name confusion,
   we introduce the following notation which should be consistent with the use
   in assembly loops

   source:  source of adjacency link
   target:  target of adjacency link
 
   E.g. the cell-node adjacency for FEM assembly links  a number of
   cells with a collection of nodes.  The cells are the surces,
   and the targets are the nodes. 

   getindex(a,i,isource) aka a[i,isource]: return i-th target of  source j
   nsources(a): overall number of sources, e.g. number of cells
   ntargets(a): overall number of targets
   ntargets(a,isource): number of targets for source given by isource
   nlinks(a): number of links aka nonzero entries of adjacency matrix
   show(a): print stuff
"""

struct VariableColumnAdjacency{T}
    colentries::Vector{T}
    colstart::Vector{T}
end

VariableColumnAdjacency(t::Type{T}) where T=VariableColumnAdjacency{T}(Vector{T}(undef,0),[one(T)])
VariableColumnAdjacency()=VariableColumnAdjacency(Int64)
VariableColumnAdjacency(m::Matrix{T}) where T=VariableColumnAdjacency{T}(vec(m),collect(1:size(m,1):size(m,1)*size(m,2)+1))

function Base.show(adj::VariableColumnAdjacency)
    for isource=1:nsources(adj)
        for itarget=1:ntargets(adj,isource)
            print(adj[itarget,isource], " ")
        end
        println()
    end
end

Base.getindex(adj::VariableColumnAdjacency,i,isource)=adj.colentries[adj.colstart[isource]+i-1]
ntargets(adj::VariableColumnAdjacency,isource)=adj.colstart[isource+1]-adj.colstart[isource]
nsources(adj::VariableColumnAdjacency)=length(adj.colstart)-1
ntargets(adj::VariableColumnAdjacency)=maximum(adj.colentries)
nlinks(adj::VariableColumnAdjacency)=length(adj.colentries)

function append!(adj::VariableColumnAdjacency,column)
    for i=1:length(column)
        push!(adj.colentries,column[i])
    end
    push!(adj.colstart,length(adj.colentries)+1)
end


ntargets(adj::Matrix,isource)=size(adj)[1]
nsources(adj::Matrix)=size(adj)[2]
ntargets(adj::Matrix)=maximum(vec(adj))
nlinks(adj::Matrix)=length(adj)

const Adjacency{T}=Union{Matrix{T},VariableColumnAdjacency{T}}


function atranspose(adj::Adjacency{T}) where T
    # 0th pass: calculate number of rows !!! todo: how to call ?
    t_adj=VariableColumnAdjacency(zeros(T,nlinks(adj)),zeros(T,ntargets(adj)+1))
    
    # 1st pass: calculate new column sizes, store them in t_adj.colstart
    for isource=1:nsources(adj)
        for itarget=1:ntargets(adj,isource)
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
    for isource=2:nsources(t_adj)
	save=t_adj.colstart[isource];
	t_adj.colstart[isource]=t_adj.colstart[isource-1]+delta; 
	if delta>0
            t_adj.colentries[t_adj.colstart[isource]-1]=t_adj.colstart[isource-1];
        end
	delta=save;
    end
    isource=nsources(t_adj)+1
    t_adj.colstart[isource]=t_adj.colstart[isource-1]+delta;
    if delta>0
        t_adj.colentries[t_adj.colstart[isource]-1]=t_adj.colstart[isource-1];
    end
    
    # 3rd pass: assemble new columns
    for isource=1:nsources(adj)
        for itarget=1:ntargets(adj,isource)
	    asm_idx=t_adj.colstart[adj[itarget,isource]+1]-1;
	    asm_loc=t_adj.colentries[asm_idx];
	    t_adj.colentries[asm_idx]+=1; 
	    t_adj.colentries[asm_loc]=isource;
        end
    end
    return t_adj
end

