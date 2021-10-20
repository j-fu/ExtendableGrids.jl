"""
````
function simplexgrid(coord::Array{Tc,2},
                     cellnodes::Array{Ti,2},
                     cellregions::Array{Ti,1},
                     bfacenodes::Array{Ti,2},
                     bfaceregions::Array{Ti,1}
                     ) where {Tc,Ti}
````

    Create simplex grid from five arrays.
"""
function simplexgrid(coord::Array{Tc,2},
                     cellnodes::Array{Ti,2},
                     cellregions::Array{Ti,1},
                     bfacenodes::Array{Ti,2},
                     bfaceregions::Array{Ti,1}
                     ) where {Tc,Ti}
    @assert size(coord,2)>0
    dim=size(coord,1)
    if dim==1
        eltype=Edge1D
        btype=Vertex0D
        csys=Cartesian1D
    elseif dim==2
        eltype=Triangle2D
        btype=Edge1D
        csys=Cartesian2D
    elseif dim==3
        eltype=Tetrahedron3D
        btype=Triangle2D
        csys=Cartesian3D
    end
    
    grid=ExtendableGrid{Tc,Ti}()
    grid[Coordinates]=coord
    grid[CellNodes]=cellnodes
    grid[CellRegions]=cellregions
    grid[CellGeometries]=VectorOfConstants(eltype,length(cellregions))
    grid[BFaceNodes]=bfacenodes
    grid[BFaceRegions]=bfaceregions
    grid[BFaceGeometries]=VectorOfConstants(btype,length(bfaceregions))
    grid[CoordinateSystem]=csys
    return grid
end


"""
````
function simplexgrid(coord::Array{Tc,2},
                     cellnodes::Array{Ti,2},
                     cellregions::Array{Ti,1},
                     bfacenodes::Array{Ti,2},
                     bfaceregions::Array{Ti,1}
                     bedgenodes::Array{Ti,2},
                     bedgeregions::Array{Ti,1}
                     ) where {Tc,Ti}
````

    Create simplex grid from coordinates, cell-nodes-adjancency, cell-region-numbers,
    boundary-face-nodes adjacency, boundary-face-region-numbers, boundary-edge-nodes, and
    boundary-edge-region-numbers arrays.
"""
function simplexgrid(coord::Array{Tc,2},
                     cellnodes::Array{Ti,2},
                     cellregions::Array{Ti,1},
                     bfacenodes::Array{Ti,2},
                     bfaceregions::Array{Ti,1},
                     bedgenodes::Array{Ti,2},
                     bedgeregions::Array{Ti,1}
                     ) where {Tc,Ti}
    grid = simplexgrid(coord, cellnodes, cellregions, bfacenodes, bfaceregions)
    grid[BEdgeNodes]   = bedgenodes
    grid[BEdgeRegions] = bedgeregions
    return grid
end


simplexgrid(C,CN,CR,BFN,BFR)=simplexgrid(collect(C),collect(CN),collect(CR),collect(BFN), collect(BFR))


##########################################################
abstract type XCoordinates <: AbstractGridFloatArray1D end
abstract type YCoordinates <: AbstractGridFloatArray1D end
abstract type ZCoordinates <: AbstractGridFloatArray1D end


"""
$(TYPEDSIGNATURES)

Constructor for 1D grid.

Construct 1D grid from an array of node cordinates.
It creates two boundary regions with index 1 at the left end and
index 2 at the right end.

Primal grid holding unknowns: marked by `o`, dual
grid marking control volumes: marked by `|`.

```@raw html
 o-----o-----o-----o-----o-----o-----o-----o-----o
 |--|-----|-----|-----|-----|-----|-----|-----|--|
```
"""
function simplexgrid(_X::AbstractVector)
    X=collect_or_assign(_X)
    #    is_monotone(X) || error("X not monotone")
    coord=reshape(X,1,length(X))
    cellnodes=zeros(Int32,2,length(X)-1)
    cellregions=zeros(Int32,length(X)-1)
    for i=1:length(X)-1 
        cellnodes[1,i]=i
        cellnodes[2,i]=i+1
        cellregions[i]=1
    end
    bfacenodes=Array{Int32}(undef,1,2)
    bfaceregions=zeros(Int32,2)
    bfacenodes[1,1]=1
    bfacenodes[1,2]=length(X)
    bfaceregions[1]=1
    bfaceregions[2]=2
    grid=simplexgrid(coord,
                     cellnodes,
                     cellregions,
                     bfacenodes,
                     bfaceregions)
    grid[XCoordinates]=X
    grid
end


##########################################################
"""
$(TYPEDSIGNATURES)

Constructor for 2D grid
from coordinate arrays. 
Boundary region numbers count counterclockwise:

| location  |  number |
| --------- | ------- |
| south     |       1 |
| east      |       2 |
| north     |       3 |
| west      |       4 |

"""
function  simplexgrid(_X::AbstractVector, _Y::AbstractVector)
    X=collect_or_assign(_X)
    Y=collect_or_assign(_Y)
    is_monotone(X) || error("X not monotone")
    is_monotone(Y) || error("Y not monotone")
    
    
    Tc=promote_type(eltype(X),eltype(Y))

    function leq(x, x1, x2)
        if (x>x1)
            return false
        end
        if (x>x2)
            return false
        end
        return true
    end
    
    function geq(x, x1, x2)
        if (x<x1)
            return false
        end
        if (x<x2)
            return false
        end
        return true
    end

    nx=length(X)
    ny=length(Y)
    
    hmin=X[2]-X[1]
    for i=1:nx-1
        h=X[i+1]-X[i]
        if h <hmin
            hmin=h
        end
    end
    for i=1:ny-1
        h=Y[i+1]-Y[i]
        if h <hmin
            hmin=h
        end
    end
    
    @assert(hmin>0.0)
    eps=1.0e-5*hmin

    x1=X[1]+eps
    xn=X[nx]-eps
    y1=Y[1]+eps
    yn=Y[ny]-eps
    
    
    function  check_insert_bface(ibface,coord,n1,n2)
        if (geq(x1,coord[1,n1],coord[1,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=4
        elseif (leq(xn,coord[1,n1],coord[1,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=2
        elseif (geq(y1,coord[2,n1],coord[2,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=1
        elseif (leq(yn,coord[2,n1],coord[2,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=3
        end
        ibface
    end
    
    
    num_nodes=nx*ny
    num_cells=2*(nx-1)*(ny-1)
    num_bfacenodes=2*(nx-1)+2*(ny-1)
    
    coord=zeros(Tc,2,num_nodes)
    cellnodes=zeros(Int32,3,num_cells)
    cellregions=zeros(Int32,num_cells)
    bfacenodes=Array{Int32}(undef,2,num_bfacenodes)
    bfaceregions=zeros(Int32,num_bfacenodes)
    
    ipoint=0
    for iy=1:ny
        for ix=1:nx
            ipoint=ipoint+1
            coord[1,ipoint]=X[ix]
            coord[2,ipoint]=Y[iy]
        end
    end
    @assert(ipoint==num_nodes)
    
    icell=0
    for iy=1:ny-1
        for ix=1:nx-1
	    ip=ix+(iy-1)*nx
	    p00 = ip
	    p10 = ip+1
	    p01 = ip  +nx
	    p11 = ip+1+nx
            
            icell=icell+1
            cellnodes[1,icell]=p00
            cellnodes[2,icell]=p10
            cellnodes[3,icell]=p11
            cellregions[icell]=1
            
            
            icell=icell+1
            cellnodes[1,icell]=p11
            cellnodes[2,icell]=p01
            cellnodes[3,icell]=p00
            cellregions[icell]=1
        end
    end
    @assert(icell==num_cells)
    
    #lazy way to  create boundary grid

    ibface=0
    for icell=1:num_cells
        n1=cellnodes[1,icell]
	n2=cellnodes[2,icell]
	n3=cellnodes[3,icell]
        ibface=check_insert_bface(ibface,coord,n1,n2)
	ibface=check_insert_bface(ibface,coord,n1,n3)
	ibface=check_insert_bface(ibface,coord,n2,n3)
    end
    @assert(ibface==num_bfacenodes)


    grid=simplexgrid(coord,
                     cellnodes,
                     cellregions,
                     bfacenodes,
                     bfaceregions)

    grid[XCoordinates]=X
    grid[YCoordinates]=Y
    grid
end


##########################################################
"""
$(TYPEDSIGNATURES)

Constructor for 3D grid
from coordinate arrays. 
Boundary region numbers:

| location  |  number |
| --------- | ------- |
| south     |       1 |
| east      |       2 |
| north     |       3 |
| west      |       4 |
| bottom    |       5 |
| top       |       6 |

"""
function  simplexgrid(_X::AbstractVector,_Y::AbstractVector,_Z::AbstractVector)
    X=collect_or_assign(_X)
    Y=collect_or_assign(_Y)
    Z=collect_or_assign(_Z)

    is_monotone(X) || error("X not monotone")
    is_monotone(Y) || error("Y not monotone")
    is_monotone(Z) || error("Z not monotone")

    Tc=promote_type(eltype(X),eltype(Y),eltype(Z))

    
    leq(x, x1, x2, x3)=x≤x1 && x≤x2 && x≤x3
    geq(x, x1, x2, x3)=x≥x1 && x≥x2 && x≥x3

    nx=length(X)
    ny=length(Y)
    nz=length(Z)
    
    hmin=X[2]-X[1]
    for i=1:nx-1
        h=X[i+1]-X[i]
        if h <hmin
            hmin=h
        end
    end
    for i=1:ny-1
        h=Y[i+1]-Y[i]
        if h <hmin
            hmin=h
        end
    end
    
    for i=1:nz-1
        h=Z[i+1]-Z[i]
        if h <hmin
            hmin=h
        end
    end
    
    @assert(hmin>0.0)
    eps=1.0e-5*hmin
    
    x1=X[1]+eps
    xn=X[end]-eps
    
    y1=Y[1]+eps
    yn=Y[end]-eps
    
    z1=Z[1]+eps
    zn=Z[end]-eps

    
    function  check_insert_bface(ibface,coord,bfacenodes,n1,n2,n3)
        if geq(x1,coord[1,n1],coord[1,n2],coord[1,n3])
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
            bfacenodes[3,ibface]=n3
	    bfaceregions[ibface]=4
        elseif leq(xn,coord[1,n1],coord[1,n2],coord[1,n3])
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
            bfacenodes[3,ibface]=n3
	    bfaceregions[ibface]=2
        elseif geq(y1,coord[2,n1],coord[2,n2],coord[2,n3])
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
            bfacenodes[3,ibface]=n3
	    bfaceregions[ibface]=1
        elseif leq(yn,coord[2,n1],coord[2,n2],coord[2,n3])
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
            bfacenodes[3,ibface]=n3
	    bfaceregions[ibface]=3
        elseif geq(z1,coord[3,n1],coord[3,n2],coord[3,n3])
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
            bfacenodes[3,ibface]=n3
	    bfaceregions[ibface]=5
        elseif leq(zn,coord[3,n1],coord[3,n2],coord[3,n3])
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
            bfacenodes[3,ibface]=n3
	    bfaceregions[ibface]=6
        end
        ibface
    end
    
    num_nodes=nx*ny*nz
    num_cells=6*(nx-1)*(ny-1)*(nz-1)
    num_bfacenodes=4*(nx-1)*(ny-1)+4*(nx-1)*(nz-1)+4*(ny-1)*(nz-1)

    Ti=Int32
    coord=zeros(Tc,3,num_nodes)
    cellnodes=zeros(Ti,4,num_cells)
    cellregions=zeros(Ti,num_cells)
    bfacenodes=zeros(Ti,3,num_bfacenodes)
    bfaceregions=zeros(Ti,num_bfacenodes)
    
    ipoint=0
    for iz=1:nz
        for iy=1:ny
            for ix=1:nx
                ipoint=ipoint+1
                coord[1,ipoint]=X[ix]
                coord[2,ipoint]=Y[iy]
                coord[3,ipoint]=Z[iz]
            end
        end
    end

    @assert(ipoint==num_nodes)
    
    icell=0
    nxy=nx*ny
    for iz=1:nz-1
        for iy=1:ny-1
            for ix=1:nx-1
                
	        ip=ix+(iy-1)*nx+(iz-1)*nxy;
                
	        p000 = ip ;
	        p100 = ip+1;
	        p010 = ip  +nx;
	        p110 = ip+1+nx;
	        p001 = ip	+nxy;
	        p101 = ip+1	+nxy;
	        p011 = ip  +nx+nxy;
	        p111 = ip+1+nx+nxy;

                icell=icell+1; cellregions[icell]=1; @. cellnodes[:,icell]=(p000,p100,p110,p111)
                icell=icell+1; cellregions[icell]=1; @. cellnodes[:,icell]=(p000,p100,p101,p111)
                icell=icell+1; cellregions[icell]=1; @. cellnodes[:,icell]=(p000,p010,p011,p111)
                icell=icell+1; cellregions[icell]=1; @. cellnodes[:,icell]=(p000,p010,p110,p111)
                icell=icell+1; cellregions[icell]=1; @. cellnodes[:,icell]=(p000,p001,p101,p111)
                icell=icell+1; cellregions[icell]=1; @. cellnodes[:,icell]=(p000,p001,p011,p111)
            end
        end
    end
    @assert(icell==num_cells)
    #lazy way to  create boundary grid
    # lazy but easy... bad for partitioning !!!
    
    ibface=0
    for icell=1:num_cells
        n1=cellnodes[1,icell]
	n2=cellnodes[2,icell]
	n3=cellnodes[3,icell]
	n4=cellnodes[4,icell]
        ibface=check_insert_bface(ibface,coord,bfacenodes,n1,n2,n3)
	ibface=check_insert_bface(ibface,coord,bfacenodes,n1,n2,n4)
	ibface=check_insert_bface(ibface,coord,bfacenodes,n1,n3,n4)
	ibface=check_insert_bface(ibface,coord,bfacenodes,n2,n3,n4)
    end
    @assert(ibface==num_bfacenodes)
    
    
    grid=simplexgrid(coord,
                     cellnodes,
                     cellregions,
                     bfacenodes,
                     bfaceregions)

    grid[XCoordinates]=X
    grid[YCoordinates]=Y
    grid[ZCoordinates]=Z
    grid
end




######################################################
"""
$(TYPEDSIGNATURES)
  
Read grid from file. Currently for pdelib sg format only
"""
function simplexgrid(file::String;format="")
    (fbase,fext)=splitext(file)
    if format==""
        format=fext[2:end]
    end
    @assert format=="sg"
    
    tks=TokenStream(file)
    expecttoken(tks,"SimplexGrid")
    version=parse(Float64,gettoken(tks))
    version20=false;

    if (version==2.0)
        version20=true;
    elseif (version==2.1)
        version20=false;
    else
        error("Read grid: wrong format version: $(version)")
    end

    dim::Int32=0
    coord=Array{Float64,2}(undef,0,0)
    cells=Array{Int32,2}(undef,0,0)
    regions=Array{Int32,1}(undef,0)
    faces=Array{Int32,2}(undef,0,0)
    bregions=Array{Int32,1}(undef,0)
    while(true)
        if (trytoken(tks,"DIMENSION"))
            dim=parse(Int32,gettoken(tks));
        elseif (trytoken(tks,"NODES")) 
            nnodes=parse(Int32,gettoken(tks));
            embdim=parse(Int32,gettoken(tks));
            if(dim!=embdim)
                error("Dimension error (DIMENSION $(dim)) in section NODES")
            end
            coord=Array{Float64,2}(undef,dim,nnodes)
            for inode=1:nnodes
                for idim=1:embdim
                    coord[idim,inode]=parse(Float64,gettoken(tks))
                end
            end
        elseif (trytoken(tks,"CELLS"))
            ncells=parse(Int32,gettoken(tks));
            cells=Array{Int32,2}(undef,dim+1,ncells)
            regions=Array{Int32,1}(undef,ncells)
            for icell=1:ncells
                for inode=1:dim+1
                    cells[inode,icell]=parse(Int32,gettoken(tks));
                end
                regions[icell]=parse(Int32,gettoken(tks));
	        if version20
		    for j=1:dim+1
		        gettoken(tks);  # skip file format garbage
                    end
                end
            end
        elseif (trytoken(tks,"FACES"))
            nfaces=parse(Int32,gettoken(tks));
            faces=Array{Int32,2}(undef,dim,nfaces)
            bregions=Array{Int32,1}(undef,nfaces)
            for iface=1:nfaces
                for inode=1:dim
                    faces[inode,iface]=parse(Int32,gettoken(tks));
                end
                bregions[iface]=parse(Int32,gettoken(tks));
	        if (version20)
		    for j=1:dim+2
		        gettoken(tks); #skip file format garbage
                    end
                end
            end
        else
            expecttoken(tks,"END")
            break
        end
    end
    simplexgrid(coord,cells,regions,faces,bregions);
end

"""
$(TYPEDSIGNATURES)
  
Write grid to file. Currently for pdelib sg format only
"""
function Base.write(fname::String, g::ExtendableGrid; format="")
    (fbase,fext)=splitext(fname)
    if format==""
        format=fext[2:end]
    end
    @assert format=="sg"

    dim_g=dim_grid(g)
    dim_s=dim_space(g)
    nn=num_nodes(g)
    nc=num_cells(g)
    nbf=num_bfaces(g)
    coord=g[Coordinates]
    cellnodes=g[CellNodes]
    bfacenodes=g[BFaceNodes]
    cellregions=g[CellRegions]
    bfaceregions=g[BFaceRegions]

    # TODO: replace @sprintf by someting non-allocating
    open(fname, "w") do file
        write(file,@sprintf("SimplexGrid"))
        write(file,@sprintf(" "))
        write(file,@sprintf("2.1\n"))
        write(file,@sprintf("#created by ExtendableGrids.jl (c) J.Fuhrmann et al\n"))
        write(file,@sprintf("#mailto:{fuhrmann|streckenbach}@wias-berlin.de\n"))
        write(file,@sprintf("#%s\n",Dates.format(Dates.now(),"yyyy-mm-ddTHH-mm-SS")))
              
        write(file,@sprintf("DIMENSION\n%d\n",dim_g))
        write(file,@sprintf("NODES\n%d %d\n",nn,dim_s))
        
        for inode=1:nn
            for idim=1:dim_s
	        write(file,@sprintf("%.20e ",coord[idim,inode]))
                write(file,@sprintf("\n"))
            end
        end

        write(file,@sprintf("CELLS\n%d\n",nc))
        for icell=1:nc
            for inode=1:dim_g+1
	        write(file,@sprintf("%d ", cellnodes[inode,icell]))
            end
            write(file,@sprintf("%d\n",cellregions[icell]))
        end
        
        write(file,@sprintf("FACES\n%d\n",nbf))
        for ibface=1:nbf
            for inode=1:dim_g
	        write(file,@sprintf("%d ", bfacenodes[inode,ibface]))
            end
            write(file,@sprintf("%d\n",bfaceregions[ibface]))
        end
        write(file,@sprintf("END\n"))
        flush(file)
        flush(file)
    end
    nothing
end



"""
$(TYPEDSIGNATURES)

Merge two grids along their common boundary facets. 

- g1: First grid to be merged
- g2: Second grid to be merged
- breg:  Interior boundary region number of merged facets. If zero (default value), no extra facets are generated.
- tol:  Distance below which two points are seen as identical. Default: 1.0e-10

"""
function glue(g1,g2; breg=0, tol=1.0e-10)
    
    dim=dim_space(g1)

    bfn1=g1[BFaceNodes]
    bfn2=g2[BFaceNodes]

    nbf1=size(bfn1,2)
    nbf2=size(bfn2,2)

    coord1=g1[Coordinates]
    coord2=g2[Coordinates]

    nn1=size(coord1,2)
    nn2=size(coord2,2)

    # numbers of faces in grid1 which match nodes in grid2
    matching_faces=zeros(Int,nbf2)
    n_matching_faces=0


    # numbers of nodes in grid1 which match nodes in grid2
    matching_nodes=zeros(Int,nn2)
    n_matching_nodes=0

    
    # Add matching pair to index of matching pairs
    function add_match!(match_list, nmatch, i1,i2)
        @assert match_list[i2]==0 || match_list[i2]==i1
        if match_list[i2]==0
            match_list[i2]=i1;
            nmatch+=1
        end
        nmatch
    end

    # Add two point indices to list of matching points
    add_matching_points(ip1,ip2) = n_matching_nodes = add_match!(matching_nodes,n_matching_nodes,ip1,ip2)

    # Add two face indices to list of matching faces
    add_matching_faces(if1,if2) = n_matching_faces = add_match!(matching_faces,n_matching_faces,if1,if2)

    ## Check if two points match (distance < tol)
    points_match(ip1, ip2) = @views norm(coord1[:,ip1]-coord2[:,ip2])<tol

    # Check if points in faces match, if so, add them to point matching index
    # Check if faces match, if so, add them to face matching index
    function match_faces(if1,if2) 
        nmatch=0
        for ip1=1:dim
            for ip2=1:dim
                xip1=bfn1[ip1,if1]
                xip2=bfn2[ip2,if2]
                if points_match(xip1,xip2)
                    add_matching_points(xip1,xip2)
                    nmatch+=1
                end
            end
        end
        
      if nmatch==dim
          add_matching_faces(if1,if2);
      end
    end

    
    # Run over all pairs of boundary faces and try to match them
    for ibf1=1:nbf1
        for ibf2=1:nbf2
            match_faces(ibf1,ibf2)
        end
    end
    @info "glue: matches found: nodes $(n_matching_nodes) bfaces: $(n_matching_faces)"



    
    creg1=g1[CellRegions]
    creg2=g2[CellRegions]
    
    breg1=g1[BFaceRegions]
    breg2=g2[BFaceRegions]
    
    cn1=g1[CellNodes]
    cn2=g2[CellNodes]

    nc1=size(cn1,2)
    nc2=size(cn2,2)


    # transposed list of matching faces
    mf1=zeros(Int,nbf1)
    for if2=1:nbf2
        if matching_faces[if2]!=0
            mf1[matching_faces[if2]]=if2
        end
    end

    mfac = breg == 0 ? 2 : 1

    coordx=zeros(dim,nn1+nn2-n_matching_nodes)
    cnx=zeros(Int, dim+1,nc1+nc2)
    cregx=zeros(Int,nc1+nc2)
    bfnx=zeros(Int,dim,nbf1+nbf2-mfac*n_matching_faces)
    bregx=zeros(Int,nbf1+nbf2-mfac*n_matching_faces)



    #   copy all data from g1 into new fields
    for ix=1:nn1
        @views coordx[:,ix].=coord1[:,ix]
    end

    for ix=1:nc1
        cregx[ix]=creg1[ix]
        @views cnx[:,ix].=cn1[:,ix]
    end

    ibfx=1
    for ibf1=1:nbf1
        if mf1[ibf1]!=0 && breg == 0
            continue
        end
        if mf1[ibf1]!=0
            bregx[ibfx]=breg
        else
            bregx[ibfx]=breg1[ibf1]
        end
        @views bfnx[:,ibfx].=bfn1[:,ibf1]
        ibfx+=1
    end


    # re-calculate matching nodes with new global numbers
    # add missing coordinates
    ix=nn1+1
    for in2=1:nn2
        if matching_nodes[in2]!=0
            continue
        end
        matching_nodes[in2]=ix
        @views coordx[:,ix].=coord2[:,in2]
        ix=ix+1
    end


    # copy missing data from g2 into arrays.
    ix=nc1+1
    for ic2=1:nc2
        cregx[ix]=creg2[ic2]
        for id=1:dim+1
            cnx[id,ix]=matching_nodes[cn2[id,ic2]]
        end
        ix=ix+1
    end

    for ibf2=1:nbf2
        if matching_faces[ibf2]!=0
            continue
        end
        bregx[ibfx]=breg2[ibf2]
        for id=1:dim
            bfnx[id,ibfx]=matching_nodes[bfn2[id,ibf2]]
        end
        ibfx+=1
    end
    @assert ibfx == nbf1+nbf2-mfac*n_matching_faces+1
    simplexgrid(coordx,cnx,cregx,bfnx,bregx)
end
