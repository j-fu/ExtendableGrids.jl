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
$(TYPEDEF)

Adjacency describing cells per boundary or interior face
"""
abstract type BFaceCells <: AbstractGridAdjacency end

"""
$(TYPEDEF)

Adjacency describing outer normals to boundary faces
"""
abstract type BFaceNormals <: AbstractGridComponent end



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
    cellnode_adj=asparse(cellnodes)
    
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
    celledge_adj=asparse(celledges)
    
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

function prepare_bfacecells!(grid)
    cn=grid[CellNodes]
    bn=grid[BFaceNodes]
    dim=dim_space(grid)
    abc=asparse(cn)'*asparse(bn)
    abcx=dropzeros!(SparseMatrixCSC(abc.m,abc.n, abc.colptr, abc.rowval, map( i-> i==dim,  abc.nzval)))
    grid[BFaceCells]=VariableTargetAdjacency(abcx)
    true
end

ExtendableGrids.instantiate(grid, ::Type{BFaceCells})=prepare_bfacecells!(grid) && grid[BFaceCells]

normal!(normal,nodes,coord,::Type{Val{1}}) = normal[1]=1.0

function normal!(normal,nodes,coord,::Type{Val{2}})
    normal[1] = -(coord[2,nodes[1]]-coord[2,nodes[2]])
    normal[2] =   coord[1,nodes[1]]-coord[1,nodes[2]]
    d=norm(normal)
    normal[1]/=d
    normal[2]/=d
end

function normal!(normal,nodes,coord,::Type{Val{3}})
    ax=coord[1,nodes[1]]-coord[1,nodes[2]]
    ay=coord[2,nodes[1]]-coord[2,nodes[2]]
    az=coord[3,nodes[1]]-coord[3,nodes[2]]
    
    bx=coord[1,nodes[1]]-coord[1,nodes[3]]
    by=coord[2,nodes[1]]-coord[2,nodes[3]]
    bz=coord[3,nodes[1]]-coord[3,nodes[3]]
    
    
    normal[1]=(ay*bz-by*az)
    normal[2]=(az*bx-bz*ax)
    normal[3]=(ax*by-bx*ay)
    
    d=norm(normal)
    normal[1]/=d
    normal[2]/=d
    normal[3]/=d
end

function midpoint!(mid,nodes,coord)
    dim=size(coord,1)
    nn=size(nodes,1)
    for i=1:dim
        mid[i]=sum(coord[i,nodes])/nn
    end
end

function adjust!(normal,cmid,bmid)
    d=dot(normal,bmid-cmid)
    if d<0.0
        normal.*= -1
    end
end

function prepare_bfacenormals!(grid)
    bfnodes=grid[BFaceNodes]
    nbf=size(bfnodes,2)
    cellnodes=grid[CellNodes]
    dim=dim_space(grid)
    bfcells=grid[BFaceCells]
    bfnormals=zeros(dim,nbf)
    coord=grid[Coordinates]
    cmid=zeros(dim)
    bmid=zeros(dim)
    for ibf=1:nbf
        icell=bfcells[1,ibf]
        @views normal!(bfnormals[:,ibf],bfnodes[:,ibf],coord,Val{dim})
        @views midpoint!(cmid,cellnodes[:,icell],coord)
        @views midpoint!(bmid,bfnodes[:,ibf],coord)
        @views adjust!(bfnormals[:,ibf],cmid,bmid)
    end
    grid[BFaceNormals]=bfnormals
    true
end

ExtendableGrids.instantiate(grid, ::Type{BFaceNormals})=prepare_bfacenormals!(grid) && grid[BFaceNormals]

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

