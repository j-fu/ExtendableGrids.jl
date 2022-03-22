using WriteVTK

# conversion from AbstractElementGeometry to WriteVTK.VTKCellTypes
VTKCellType(::Type{<:AbstractElementGeometry1D}) = VTKCellTypes.VTK_LINE
VTKCellType(::Type{<:Triangle2D})                = VTKCellTypes.VTK_TRIANGLE
VTKCellType(::Type{<:Quadrilateral2D})           = VTKCellTypes.VTK_QUAD
VTKCellType(::Type{<:Tetrahedron3D})             = VTKCellTypes.VTK_TETRA
VTKCellType(::Type{<:Hexahedron3D})              = VTKCellTypes.VTK_HEXAHEDRON

"""
$(TYPEDSIGNATURES)

exports grid and optional provided data as a vtk file

- `filename`: filename of the exported file
- `grid`: grid

Each '(key, value)' pair adds another data entry to the vtk file via WriteVTK functionality.
"""

function writeVTK(filename::String, grid::ExtendableGrid{Tc,Ti}; kwargs...) where {Tc,Ti}
    
    ncells  = num_cells(grid)   # get number of cells in grid
    coords  = grid[Coordinates] # get coordinates 
    cells   = grid[CellNodes]   # get cell-node list
    geo_dim = size(coords, 1)

    cell_geo     = grid[CellGeometries]
    cell_regions = grid[CellRegions]

    vtk_cells = Array{MeshCell,1}(undef, ncells)

    for icell in 1:ncells
        vtk_cells[icell] = MeshCell(VTKCellType(cell_geo[icell]), view(cells, :, icell))
    end 

    vtk_grid(filename, coords, vtk_cells) do vtk
        for (key, value) in kwargs
            vtk[String(key)] = value
        end
    end

end
