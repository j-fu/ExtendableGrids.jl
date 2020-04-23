function simplexgrid(coord::Array{Tc,2},
                     cellnodes::Array{Ti,2},
                     cellregions::Array{Ti,1},
                     bfacenodes::Array{Ti,2},
                     bfaceregions::Array{Ti,1}
                     ) where {Tc,Ti}
    
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
    grid[CellTypes]=VectorOfConstants(eltype,length(cellregions))
    grid[BFaceNodes]=bfacenodes
    grid[BFaceRegions]=bfaceregions
    grid[BFaceTypes]=VectorOfConstants(btype,length(bfaceregions))
    grid[CoordinateSystem]=csys
    return grid
end


"""
$(SIGNATURES)

(Try to) create a subdivision of interval (a,b) stored in the 
returned array X such that 
  - `X[1]==a, X[end]==b`
  - `(X[2]-X[1])<=ha+tol*(b-a)`
  - `(X[end]-X[end-1])<=hb+tol*(b-a)`
  - There is a number q such that  `X[i+1]-X[i] == q*(X[i]-X[i-1])`
  - X is the array with the minimal possible number of points with the above property
  
Caveat: the algorithm behind this is  well tested but unproven.

Returns an Array containing the points of the subdivision.
"""
function geomspace(a::Tv, b::Tv, ha::Tv, hb::Tv; tol=1.0e-10) where Tv
    
    function _geomspace0(l,h0, hl, tol=1.0e-10)
        
        @assert (l>0.0)
        @assert (h0>0.0)
        @assert (hl>=h0)
        @assert((hl+h0)<l)
        
        #  We need to  adjust two things:
        
        # The sum of the geometric progression must
        # match the length of the interval, so lmismatch
        # should be zero:
        function lmismatch(q,k)
            return l - h0*(1-q^k)/(1-q)
        end
        
        # The claim from experimenral evidence (Wolfram) 
        # is that, if written as a polynomial,
        # it has two real zeros: one and the value searched for which
        # is slightly larger than one. All other zeros are on one circle.
        
        # The size of the last interval should be close to 
        # to hl, so hmismatch should be close to one and not larger than one
        function  hmismatch(q,k)
            return  h0*q^(k-1)/hl
        end

        # define initial number of intervals from
        # average of minmal and maximal h
        n=Int32(ceil((2.0*l/(h0+hl))))
        
        if n==1
            n=2
        end
        
        # define initial q such that hmismatch is one.
        q=(hl/h0)^(1.0/(n-1.0))
        
        # Iteration until both mismatches are satisfactory
        # Outer loop runs until hmismatch is less than 1
        hmiss=10.0 # some initial value >1 just to run the loop at least once
        if abs(q-1.0)<tol
            hmiss=1.0
        end
        while  hmiss>1.0
            # increase number of intervals until
            # lmismatch becomes less than zero
            while  lmismatch(q,n)>0.0
                n+=1
            end
            
            # find initial interval for q containing
            # value with zero lmismatch 
            ns=0
            nsmax=1000
            
            while lmismatch(q,n)<0.0 &&  ns<nsmax
                q*=0.99
                ns+=1
            end
            
            xl=q
            xr=q/0.99
            @assert ns<nsmax
            
            # bisection to define q with zero lmismatch
            ns=0
            xm=0.5*(xl+xr)
            while (xr-xl)>tol && ns<nsmax
                ns+=1
                mmm=lmismatch(xm,n)
                if mmm==0.0
                    break
                elseif   lmismatch(xl,n)*mmm<0.0
                    xr=xm
                else
                    xl=xm
                end
                xm=0.5*(xl+xr)
            end
            q=xm
            @assert ns<nsmax
            hmiss=hmismatch(q,n)
            if hmiss>1.0 
                n=n+1
            end
        end
        #  printf("%d %g %g %g\n",n,q,lmismatch(q,n),hmismatch(q,n))
        
        X = Array{Tv,1}(undef,n+1)
        X[1]=0
        h=h0
        for i=1:n
            X[i+1]=X[i]+h
            h*=q
        end
        X[n+1]=l
        return X
    end

    # Map things to [0,b-a]
    @assert (ha>0.0)
    @assert (hb>0.0)
    @assert (a<b)


    tol=tol*(b-a)
    if ha<=hb
        X=_geomspace0(b-a,ha,hb,tol)
        X.+=a
    else
        X=-reverse(_geomspace0(b-a,hb,ha,tol))
        X.+=b
    end

    @assert (X[2]-X[1])<=ha+tol
    @assert (X[end]-X[end-1])<=hb+tol
    
    return X
end

"""
$(SIGNATURES)

Glue together two vectors a and b resulting in a vector c. They last element 
of a shall be equal (up to tol) to the first element of b.
The result fulfills `length(c)=length(a)+length(b)-1`
"""
function glue(a::Vector{Tv}, b::Vector{Tv}; tol=1.0e-10) where Tv
    #assert(is_monotone(a));
    #assert(is_monotone(b));
    na=length(a)
    nb=length(b)
    
    d=b[1]-a[na-1]
    @assert(d>0)
    d=b[1]-a[na]
    @assert(d>-tol)
    @assert(d<tol)

    c=Vector{Tv}(undef,na+nb-1)
    ic=0
    for ia=1:na
        ic+=1
        c[ic]=a[ia]
    end
    for ib=2:nb
        ic+=1
        c[ic]=b[ib]
    end
    return c
end


##########################################################
"""
$(SIGNATURES)

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
function simplexgrid(X::AbstractArray{Tc,1}) where {Tc}
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
    return simplexgrid(coord,
                cellnodes,
                cellregions,
                bfacenodes,
                bfaceregions)
end


##########################################################
"""
$(SIGNATURES)

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
function  simplexgrid(X::AbstractArray{Tc,1},Y::AbstractArray{Tc,1}) where {Tc}

    
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
    
    
    function  check_insert_bface(n1,n2)
                
        if (geq(x1,coord[1,n1],coord[1,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=4
            return
        end
        if (leq(xn,coord[1,n1],coord[1,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=2
            return
        end
        if (geq(y1,coord[2,n1],coord[2,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=1
            return
        end
        if (leq(yn,coord[2,n1],coord[2,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=3
            return
        end
    end
    
    
    num_nodes=nx*ny
    num_cells=2*(nx-1)*(ny-1)
    num_bfacenodes=2*(nx-1)+2*(ny-1)
    
    coord=zeros(Tc,2,num_nodes)
    cellnodes=zeros(Int32,3,num_cells)
    cellregions=zeros(Int32,num_cells)
    bfacenodes=Array{Int32}(undef,2,num_bfacenodes)
#    resize!(bfacenodes,2,num_bfacenodes)
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
        check_insert_bface(n1,n2)
	check_insert_bface(n1,n3)
	check_insert_bface(n2,n3)
    end
    @assert(ibface==num_bfacenodes)


    return simplexgrid(coord,
                cellnodes,
                cellregions,
                bfacenodes,
                bfaceregions)
end




"""
Create Grid from Triangle input data.
"""
function simplexgrid(flags::String, input::Triangulate.TriangulateIO)

    triout,vorout=Triangulate.triangulate(flags,input)

    pointlist=triout.pointlist
    if eltype(pointlist)!=Float64
        pointlist=Array{Float64,2}(pointlist)
    end
    
    trianglelist=triout.trianglelist
    if eltype(trianglelist)!=Int32
        trianglelist=Array{Int32,2}(trianglelist)
    end

    cellregions=Vector{Int32}(vec(triout.triangleattributelist))
    
    segmentlist=triout.segmentlist
    if eltype(segmentlist)!=Int32
        segmentlist=Array{Int32,2}(segmentlist)
    end
    
    segmentmarkerlist=triout.segmentmarkerlist
    if eltype(segmentmarkerlist)!=Int32
        segmentmarkerlist=Array{Int32,2}(segmentmarkerlist)
    end

    grid=ExtendableGrid{Float64,Int32}()
    grid[Coordinates]=pointlist
    grid[CellRegions]=cellregions
    grid[CellTypes]=VectorOfConstants(Simplex2D,length(cellregions))
    grid[BFaceRegions]=segmentmarkerlist
    grid[BFaceTypes]=VectorOfConstants(Simplex1D,length(segmentmarkerlist))
    grid[CellNodes]=trianglelist
    grid[BFaceNodes]=segmentlist
    return grid
end

"""
Create Grid from a number of input arrays.
The 2D input arrays are transposed if necessary and converted to
the proper data types for Triangulate.

This conversion is not performed if the data types are thos
indicated in the defaults and the leading dimension of 2D arrays
corresponds to the space dimension.
"""
function simplexgrid(;flags::String="pAaqDQ",
                     points=Array{Cdouble,2}(undef,0,0),
                     bfaces=Array{Cint,2}(undef,0,0),
                     bfaceregions=Array{Cint,1}(undef,0),
                     regionpoints=Array{Cdouble,2}(undef,0,0),
                     regionnumbers=Array{Cint,1}(undef,0),
                     regionvolumes=Array{Cdouble,1}(undef,0)
                  )
    @assert ndims(points)==2
    if size(points,2)==2
        points=transpose(points)
    end
    if typeof(points)!=Array{Cdouble,2}
        points=Array{Cdouble,2}(points)
    end
    @assert(size(points,2)>2)
    
    @assert ndims(bfaces)==2
    if size(bfaces,2)==2
        bfaces=transpose(bfaces)
    end
    if typeof(bfaces)!=Array{Cint,2}
        bfaces=Array{Cint,2}(bfaces)
    end
    @assert(size(bfaces,2)>0)
    
    @assert ndims(bfaceregions)==1
    @assert size(bfaceregions,1)==size(bfaces,2)
    if typeof(bfaceregions)!=Array{Cint,1}
        bfaceregions=Array{Cint,1}(bfaceregions)
    end
    
    @assert ndims(regionpoints)==2
    if size(regionpoints,2)==2
        regionpoints=transpose(regionpoints)
    end
    if typeof(regionpoints)!=Array{Cdouble,2}
        regionpoints=Array{Cdouble,2}(regionpoints)
    end
    @assert(size(regionpoints,2)>0)
    
    @assert ndims(regionnumbers)==1
    @assert ndims(regionvolumes)==1
    @assert size(regionnumbers,1)==size(regionpoints,2)
    @assert size(regionvolumes,1)==size(regionpoints,2)
    
    nholes=0
    nregions=0
    for i=1:length(regionnumbers)
        if regionnumbers[i]==0
            nholes+=1
        else
            nregions+=1
        end
    end


    
    regionlist=Array{Cdouble,2}(undef,4,nregions)
    holelist=Array{Cdouble,2}(undef,2,nholes)
    
    ihole=1
    iregion=1
    for i=1:length(regionnumbers)
        if regionnumbers[i]==0
            holelist[1,iregion]=regionpoints[1,i]
            holeist[2,iregion]=regionpoints[2,i]
            ihole+=1
        else
            regionlist[1,iregion]=regionpoints[1,i]
            regionlist[2,iregion]=regionpoints[2,i]
            regionlist[3,iregion]=regionnumbers[i]
            regionlist[4,iregion]=regionvolumes[i]
            iregion+=1
        end
    end
    tio=Triangulate.TriangulateIO()
    tio.pointlist=points
    tio.segmentlist=bfaces
    tio.segmentmarkerlist=bfaceregions
    tio.regionlist=regionlist
    tio.holelist=holelist
    return generate(flags,tio)
end


