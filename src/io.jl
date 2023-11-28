using WriteVTK

# conversion from AbstractElementGeometry to WriteVTK.VTKCellTypes
VTKCellType(::Type{<:AbstractElementGeometry1D}) = VTKCellTypes.VTK_LINE
VTKCellType(::Type{<:Triangle2D}) = VTKCellTypes.VTK_TRIANGLE
VTKCellType(::Type{<:Quadrilateral2D}) = VTKCellTypes.VTK_QUAD
VTKCellType(::Type{<:Tetrahedron3D}) = VTKCellTypes.VTK_TETRA
VTKCellType(::Type{<:Hexahedron3D}) = VTKCellTypes.VTK_HEXAHEDRON

"""
$(TYPEDSIGNATURES)

exports grid and optional provided data as a vtk file

- `filename`: filename of the exported file
- `grid`: grid

Each '(key, value)' pair adds another data entry to the vtk file via WriteVTK functionality.
"""

function writeVTK(filename::String, grid::ExtendableGrid{Tc, Ti}; kwargs...) where {Tc, Ti}
    ncells = num_cells(grid)   # get number of cells in grid
    coords = grid[Coordinates] # get coordinates 
    cells = grid[CellNodes]   # get cell-node list
    geo_dim = size(coords, 1)

    cell_geo = grid[CellGeometries]

    vtk_cells = Array{MeshCell, 1}(undef, ncells)

    for icell = 1:ncells
        vtk_cells[icell] = MeshCell(VTKCellType(cell_geo[icell]), view(cells, :, icell))
    end

    vtk_grid(filename, coords, vtk_cells) do vtk
        for (key, value) in kwargs
            vtk[String(key)] = value
        end
    end
end

"""
$(TYPEDSIGNATURES)
  
Write grid to file. Currently for pdelib sg and Gmsh formats.
"""
function Base.write(fname::String, g::ExtendableGrid; format = "")
    (fbase, fext) = splitext(fname)
    if format == ""
        format = fext[2:end]
    end
    if format == "msh"
        try
            simplexgrid_to_gmsh(g; filename = fname)
        catch e
            throw(ErrorException("Missing Gmsh extension. Add Gmsh.jl to your environment and import it to write $(fext) files."))
        end
        return
    end
    @assert format == "sg"

    dim_g = dim_grid(g)
    dim_s = dim_space(g)
    nn = num_nodes(g)
    nc = num_cells(g)
    nbf = num_bfaces(g)
    coord = g[Coordinates]
    cellnodes = g[CellNodes]
    bfacenodes = g[BFaceNodes]
    cellregions = g[CellRegions]
    bfaceregions = g[BFaceRegions]

    # TODO: replace @sprintf by someting non-allocating
    open(fname, "w") do file
        write(file, @sprintf("SimplexGrid"))
        write(file, @sprintf(" "))
        write(file, @sprintf("2.1\n"))
        write(file, @sprintf("#created by ExtendableGrids.jl (c) J.Fuhrmann et al\n"))
        write(file, @sprintf("#mailto:{fuhrmann|streckenbach}@wias-berlin.de\n"))
        write(file, @sprintf("#%s\n", Dates.format(Dates.now(), "yyyy-mm-ddTHH-mm-SS")))

        write(file, @sprintf("DIMENSION\n%d\n", dim_g))
        write(file, @sprintf("NODES\n%d %d\n", nn, dim_s))

        for inode = 1:nn
            for idim = 1:dim_s
                write(file, @sprintf("%.20e ", coord[idim, inode]))
                write(file, @sprintf("\n"))
            end
        end

        write(file, @sprintf("CELLS\n%d\n", nc))
        for icell = 1:nc
            for inode = 1:(dim_g + 1)
                write(file, @sprintf("%d ", cellnodes[inode, icell]))
            end
            write(file, @sprintf("%d\n", cellregions[icell]))
        end

        write(file, @sprintf("FACES\n%d\n", nbf))
        for ibface = 1:nbf
            for inode = 1:dim_g
                write(file, @sprintf("%d ", bfacenodes[inode, ibface]))
            end
            write(file, @sprintf("%d\n", bfaceregions[ibface]))
        end
        write(file, @sprintf("END\n"))
        flush(file)
        flush(file)
    end
    nothing
end

######################################################
"""
$(TYPEDSIGNATURES)
  
Read grid from file. Currently for pdelib sg and Gmsh formats.
"""
function simplexgrid(file::String; format = "")
    Ti = Cint
    (fbase, fext) = splitext(file)
    if format == ""
        format = fext[2:end]
    end
    if format == "msh" || format == "geo"
        grid = nothing
        try
            grid = simplexgrid_from_gmsh(file)
        catch e
            throw(ErrorException("Missing Gmsh extension. Add Gmsh.jl to your environment and import it to read $(fext) files."))
        end
        return grid
    end

    @assert format == "sg"

    tks = TokenStream(file)
    expecttoken(tks, "SimplexGrid")
    version = parse(Float64, gettoken(tks))
    version20 = false

    if (version == 2.0)
        version20 = true
    elseif (version == 2.1)
        version20 = false
    else
        error("Read grid: wrong format version: $(version)")
    end

    dim::Ti = 0
    coord = Array{Float64, 2}(undef, 0, 0)
    cells = Array{Ti, 2}(undef, 0, 0)
    regions = Array{Ti, 1}(undef, 0)
    faces = Array{Ti, 2}(undef, 0, 0)
    bregions = Array{Ti, 1}(undef, 0)
    while (true)
        if (trytoken(tks, "DIMENSION"))
            dim = parse(Ti, gettoken(tks))
        elseif (trytoken(tks, "NODES"))
            nnodes = parse(Ti, gettoken(tks))
            embdim = parse(Ti, gettoken(tks))
            if (dim != embdim)
                error("Dimension error (DIMENSION $(dim)) in section NODES")
            end
            coord = Array{Float64, 2}(undef, dim, nnodes)
            for inode = 1:nnodes
                for idim = 1:embdim
                    coord[idim, inode] = parse(Float64, gettoken(tks))
                end
            end
        elseif (trytoken(tks, "CELLS"))
            ncells = parse(Ti, gettoken(tks))
            cells = Array{Ti, 2}(undef, dim + 1, ncells)
            regions = Array{Ti, 1}(undef, ncells)
            for icell = 1:ncells
                for inode = 1:(dim + 1)
                    cells[inode, icell] = parse(Ti, gettoken(tks))
                end
                regions[icell] = parse(Ti, gettoken(tks))
                if version20
                    for j = 1:(dim + 1)
                        gettoken(tks)  # skip file format garbage
                    end
                end
            end
        elseif (trytoken(tks, "FACES"))
            nfaces = parse(Ti, gettoken(tks))
            faces = Array{Ti, 2}(undef, dim, nfaces)
            bregions = Array{Ti, 1}(undef, nfaces)
            for iface = 1:nfaces
                for inode = 1:dim
                    faces[inode, iface] = parse(Ti, gettoken(tks))
                end
                bregions[iface] = parse(Ti, gettoken(tks))
                if (version20)
                    for j = 1:(dim + 2)
                        gettoken(tks) #skip file format garbage
                    end
                end
            end
        else
            expecttoken(tks, "END")
            break
        end
    end
    simplexgrid(coord, cells, regions, faces, bregions)
end

function simplexgrid_from_gmsh end

function simplexgrid_to_gmsh end

function mixedgrid_from_gmsh end

function mixedgrid_to_gmsh end
