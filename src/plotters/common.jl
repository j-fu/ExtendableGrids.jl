# TODO: use distinguishable colors
# http://juliagraphics.github.io/Colors.jl/stable/colormapsandcolorscales/#Generating-distinguishable-colors-1
"""
$(SIGNATURES)

Color scale for grid colors.
"""
function frgb(Plotter,i,max;pastel=false)
    x=Float64(i-1)/Float64(max)
    if (x<0.5)
        r=1.0-2.0*x
        g=2.0*x
        b=0.0
    else
        r=0.0
        g=2.0-2.0*x
        b=2.0*x-1.0
    end
    if pastel
        r=0.5+0.5*r
        g=0.5+0.5*g
        b=0.5+0.5*b
    end
    if ispyplot(Plotter)
        return (r,g,b)
    end
    if ismakie(Plotter)
        return RGB(r,g,b)
    end
    if isplots(Plotter)
        return Plotter.RGB(r,g,b)
    end
end


"""
$(SIGNATURES)

Extract visible tetrahedra - those intersecting with the planes
`x=xyzcut[1]` or `y=xyzcut[2]`  or `z=xyzcut[3]`. 

Return corresponding points and facets for each region for drawing as mesh (Makie,MeshCat)
or trisurf (pyplot)
"""
function extract_visible_cells3D(grid::ExtendableGrid,xyzcut)
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    cellregions=grid[CellRegions]
    nregions=grid[NumCellRegions]
    
    cutcoord=zeros(3*4)

    function take(coord,simplex,xyzcut)
        icut=1
        for idim=1:3
            for inode=1:4
                cutcoord[icut]=coord[idim,simplex[inode]]-xyzcut[idim]
                icut=icut+1
            end
        end
        mapreduce(a->a<=0,*,cutcoord)|| mapreduce(a->a>=0,*,cutcoord)
    end

    pmark=[zeros(Int32,size(coord,2)) for ireg=1:nregions]
    faces=[ElasticArray{Int32}(undef,3,0) for iregion=1:nregions]
    npoints=zeros(Int32,nregions)

    mtet=zeros(4)
    for itet=1:size(cellnodes,2)
        iregion=cellregions[itet]
        tet=view(cellnodes,:, itet)
        if take(coord,tet,xyzcut)
            for inode=1:4
                if pmark[iregion][tet[inode]]==0
                    npoints[iregion]+=1
                    pmark[iregion][tet[inode]]=npoints[iregion]
                end
                mtet[inode]=pmark[iregion][tet[inode]]
            end
            append!(faces[iregion],(mtet[1],mtet[2],mtet[3]))
            append!(faces[iregion],(mtet[1],mtet[2],mtet[4]))
            append!(faces[iregion],(mtet[2],mtet[3],mtet[4]))
            append!(faces[iregion],(mtet[3],mtet[1],mtet[4]))
        end
    end
    
    points=[ Array{Float32,2}(undef,3,npoints[iregion]) for iregion=1:nregions]
    for i=1:size(coord,2)
        for iregion=1:nregions
            if pmark[iregion][i]>0
                @views  points[iregion][:,pmark[iregion][i]].=(coord[1,i],coord[2,i],coord[3,i])
            end
        end
    end
    points,faces
end



function extract_visible_bfaces3D(grid::ExtendableGrid,xyzcut)
    cutcoord=zeros(3)
    
    function take(coord,simplex,xyzcut)
        for idim=1:3
            for inode=1:3
                cutcoord[inode]=coord[idim,simplex[inode]]-xyzcut[idim]
            end
            if !mapreduce(a->a<=0,*,cutcoord)
                return false
            end
        end
        return true
    end
    
    coord=grid[Coordinates]
    nbfaces=num_bfaces(grid)
    bfacenodes=grid[BFaceNodes]
    bfaceregions=grid[BFaceRegions]
    nbregions=grid[NumBFaceRegions]

    pmark=[zeros(Int32,size(coord,2)) for ireg=1:nbregions]
    faces=[ElasticArray{Int32}(undef,3,0) for iregion=1:nbregions]
    npoints=zeros(Int32,nbregions)
    

    for i=1:nbfaces
        ibregion=bfaceregions[i]
        tri=view(bfacenodes,:, i)
        if take(coord,tri,xyzcut)
            for inode=1:3
                if pmark[ibregion][tri[inode]]==0
                    npoints[ibregion]+=1
                    pmark[ibregion][tri[inode]]=npoints[ibregion]
                end
            end
            tri=map(i->pmark[ibregion][i],tri)
            append!(faces[ibregion],tri)
        end
    end
    points=[ Array{Float32,2}(undef,3,npoints[ibregion]) for ibregion=1:nbregions]
    for i=1:size(coord,2)
        for ibregion=1:nbregions
            if pmark[ibregion][i]>0
                @views points[ibregion][:,pmark[ibregion][i]].=(coord[1,i],coord[2,i],coord[3,i])
            end
        end
    end
    points,faces
end


function extract_visible_bfaces3D(grid::ExtendableGrid,func,xyzcut)
    cutcoord=zeros(3)
    
    function take(coord,simplex,xyzcut)
        for idim=1:3
            for inode=1:3
                cutcoord[inode]=coord[idim,simplex[inode]]-xyzcut[idim]
            end
            if !mapreduce(a->a<=0,*,cutcoord)
                return false
            end
        end
        return true
    end
    
    coord=grid[Coordinates]
    nbfaces=num_bfaces(grid)
    bfacenodes=grid[BFaceNodes]
    
    pmark=zeros(Int32,size(coord,2))
    faces=ElasticArray{Int32}(undef,3,0)
    npoints=0
    
    for i=1:nbfaces
        tri=view(bfacenodes,:, i)
        if take(coord,tri,xyzcut)
            for inode=1:3
                if pmark[tri[inode]]==0
                    npoints+=1
                    pmark[tri[inode]]=npoints
                end
            end
            tri=map(i->pmark[i],tri)
            append!(faces,tri)
        end
    end
    
    points=Array{Float32,2}(undef,3,npoints)
    values=Vector{Float32}(undef,npoints)
    for i=1:size(coord,2)
        if pmark[i]>0
            @views points[:,pmark[i]].=(coord[1,i],coord[2,i],coord[3,i])
            values[pmark[i]]=func[i]
        end
    end
    points,faces,values
end
