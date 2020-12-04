"""
$(TYPEDEF)

Adjacency describing edges per grid cell
"""
abstract type CellEdges  <: AbstractGridAdjacency end

"""
$(TYPEDEF)

Adjacency describing cells per grid edge
"""
abstract type EdgeCells  <: AbstractGridAdjacency end


"""
$(TYPEDEF)

Adjacency describing nodes per grid edge
"""
abstract type EdgeNodes <: AbstractGridAdjacency end



"""
$(SIGNATURES)

Prepare edge adjacencies (celledges, edgecells, edgenodes)

Currently depends on ExtendableSparse, we may want to remove this
adjacency.
""" 
function prepare_edges!(grid::ExtendableGrid)
    Ti=eltype(grid[CellNodes])
    cellnodes=grid[CellNodes]
    geom=grid[CellGeometries][1]
    # Create cell-node incidence matrix
    ext_cellnode_adj=ExtendableSparseMatrix{Ti,Ti}(num_nodes(grid),num_cells(grid))
    for icell=1:num_cells(grid)
        for inode=1:num_nodes(geom)
            ext_cellnode_adj[cellnodes[inode,icell],icell]=1
        end
    end
    flush!(ext_cellnode_adj)
    # Get SparseMatrixCSC from the ExtendableMatrix
    cellnode_adj=ext_cellnode_adj.cscmatrix
    
    # Create node-node incidence matrix for neigboring
    # nodes. 
    nodenode_adj=cellnode_adj*transpose(cellnode_adj)

    # To get unique edges, we set the lower triangular part
    # including the diagonal to 0
    for icol=1:length(nodenode_adj.colptr)-1
        for irow=nodenode_adj.colptr[icol]:nodenode_adj.colptr[icol+1]-1
            if nodenode_adj.rowval[irow]>=icol
                nodenode_adj.nzval[irow]=0
            end
        end
    end
    dropzeros!(nodenode_adj)


    # Now we know the number of edges and
    nedges=length(nodenode_adj.nzval)

    
    if dim_space(grid)==2
        # Let us do the Euler test (assuming no holes in the domain)
        v=num_nodes(grid)
        e=nedges
        f=num_cells(grid)+1
        @assert v-e+f==2
    end

    if dim_space(grid)==1
        @assert nedges==num_cells(grid)
    end
    
    # Calculate edge nodes and celledges
    edgenodes=zeros(Ti,2,nedges)
    celledges=zeros(Ti,num_edges(geom),num_cells(grid))
    cen=local_celledgenodes(geom)
    
    for icell=1:num_cells(grid)
        for iedge=1:num_edges(geom)
            n1=cellnodes[cen[1,iedge],icell]
            n2=cellnodes[cen[2,iedge],icell]

            # We need to look in nodenod_adj for upper triangular part entries
            # therefore, we need to swap accordingly before looking
	    if (n1<n2)
		n0=n1
		n1=n2
		n2=n0;
	    end
            
            for irow=nodenode_adj.colptr[n1]:nodenode_adj.colptr[n1+1]-1
                if nodenode_adj.rowval[irow]==n2
                    # If the coresponding entry has been found, set its
                    # value. Note that this introduces a different edge orientation
                    # compared to the one found locally from cell data
                    celledges[iedge,icell]=irow
                    edgenodes[1,irow]=n1
                    edgenodes[2,irow]=n2
                end
            end
        end
    end


    # Create sparse incidence matrix for the cell-edge adjacency
    ext_celledge_adj=ExtendableSparseMatrix{Ti,Ti}(nedges,num_cells(grid))
    for icell=1:num_cells(grid)
        for iedge=1:num_edges(geom)
            ext_celledge_adj[celledges[iedge,icell],icell]=1
        end
    end
    flush!(ext_celledge_adj)
    celledge_adj=ext_celledge_adj.cscmatrix

    # The edge cell matrix is the transpose
    edgecell_adj=SparseMatrixCSC(transpose(celledge_adj))

    # Get the adjaency array from the matrix
    edgecells=zeros(Ti,2,nedges) ## for 3D we need more here!
    for icol=1:length(edgecell_adj.colptr)-1
        ii=1
        for irow=edgecell_adj.colptr[icol]:edgecell_adj.colptr[icol+1]-1
            edgecells[ii,icol]=edgecell_adj.rowval[irow]
            ii+=1
        end
    end
    
    grid[EdgeCells]=edgecells
    grid[CellEdges]=celledges
    grid[EdgeNodes]=edgenodes
    true
end

ExtendableGrids.instantiate(grid, ::Type{CellEdges})=prepare_edges!(grid) && grid[CellEdges]
ExtendableGrids.instantiate(grid, ::Type{EdgeCells})=prepare_edges!(grid) && grid[EdgeCells]
ExtendableGrids.instantiate(grid, ::Type{EdgeNodes})=prepare_edges!(grid) && grid[EdgeNodes]


# The following methods are uses in VoronoiFVM.
# Do we need a more systematic approach here ?
"""
$(SIGNATURES)

Number of edges in grid.
"""
num_edges(grid::ExtendableGrid)=haskey(grid,EdgeNodes) ?  num_sources(grid[EdgeNodes]) : 0

"""
$(SIGNATURES)

Number of nodes of 0D vertex
"""
num_nodes(::Type{Vertex0D})=1

"""
$(SIGNATURES)

Number of edges of 0D vertex
"""
num_edges(::Type{Vertex0D})=0

const cen_Edge1D=reshape([1 2],:,1)
"""
$(SIGNATURES)

Cell-edege node numbering for 1D edge
"""
local_celledgenodes(::Type{Edge1D})=cen_Edge1D

"""
$(SIGNATURES)

Number of nodes for 1D edge
"""
num_nodes(::Type{Edge1D})=2

"""
$(SIGNATURES)

Number of edges of 1D edge
"""
num_edges(::Type{Edge1D})=1

const cen_Triangle2D=[ 2 3 1; 3 1 2]

"""
$(SIGNATURES)

Cell-edege node numbering for 2D triangle
"""
local_celledgenodes(::Type{Triangle2D})=cen_Triangle2D


const cen_Tetrahedron3D=[ 3 4 2  1 1 1; 4 2 3  4 3 2]

"""
$(SIGNATURES)

Cell-edege node numbering for 2D triangle
"""
local_celledgenodes(::Type{Tetrahedron3D})=cen_Tetrahedron3D


"""
$(SIGNATURES)

Number of nodes in 2D triangle
"""
num_nodes(::Type{Triangle2D})=3

"""
$(SIGNATURES)

Number of nodes in 3D tetrahedron
"""
num_nodes(::Type{Tetrahedron3D})=4

"""
$(SIGNATURES)

Number of edges in 2D triangle
"""
num_edges(::Type{Triangle2D})=3

"""
$(SIGNATURES)

Number of edges in 3D tetrahedron
"""
num_edges(::Type{Tetrahedron3D})=6

