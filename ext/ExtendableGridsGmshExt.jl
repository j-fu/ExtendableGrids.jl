module ExtendableGridsGmshExt

if isdefined(Base, :get_extension)
    import Gmsh: gmsh
else
    import ..Gmsh: gmsh
end


import ExtendableGrids: simplexgrid_from_gmsh, load_simplexgrid_to_gmsh
import ExtendableGrids: mixedgrid_from_gmsh, load_mixedgrid_to_gmsh

import ExtendableGrids: ExtendableGrid, simplexgrid, VariableTargetAdjacency, num_sources
import ExtendableGrids: Coordinates, CellNodes, CellRegions, BFaceNodes, BFaceRegions, CellGeometries, BFaceGeometries
import ExtendableGrids: Edge1D, Triangle2D, Quadrilateral2D, Tetrahedron3D, Hexahedron3D, Prism3D, VectorOfConstants, ElementGeometries, Cartesian2D, Cartesian3D, CoordinateSystem, num_cells


#!!! Make a license warning at initialization ? Gmsh is GPL - mention this in the readme.

###??? Do we really need this dependency here ? I would rather like to live without, the more that it seems to make some problems.
#using ExtendableGrids
using StatsBase: countmap
using Bijections


"""
````
function simplexgrid_from_gmsh(filename::String; incomplete=false)
````

Then the msh file is read and a SimplexGrid is created.
The mesh can also contain an incomplete grid. For this, the function has to be called with ````incomplete=true````.
'incomplete' means that the grid only consists of nodes and cells, it does not have a boundary. We also do not try to read the physical groups for those grids.


"""
function simplexgrid_from_gmsh(filename::String; incomplete=false)
	gmshfile_to_simplexgrid(filename; incomplete)
end

"""
````
function mixedgrid_from_gmsh(filename::String)
````

Then the msh file is read and an ExtendableGrid is created.
This only works for dim=2 grids and the orientation may be wrong.

"""
function mixedgrid_from_gmsh(filename::String)
	gmshfile_to_mixedgrid(filename)
end



"""
````
simplexgrid_from_gmsh(mod::Module; incomplete=false)
````

Then the mesh contained in the gmsh module is converted to a SimplexGrid.
The mesh can also contain an incomplete grid. For this, the function has to be called with ````incomplete=true````.
'incomplete' means that the grid only consists of nodes and cells, it does not have a boundary. We also do not try to read the physical groups for those grids.

"""
function simplexgrid_from_gmsh(mod::Module; incomplete=false)
    if !incomplete
	    mod_to_simplexgrid(mod)
	else
		incomplete_mod_to_simplexgrid(mod)
	end
end

"""
````
mixedgrid_from_gmsh(mod::Module)
````

Then the mesh contained in the gmsh module is converted to an ExtendableGrid.

"""
function mixedgrid_from_gmsh(mod::Module)
    mod_to_mixedgrid(mod)
end

"""
````
function load_simplexgrid_to_gmsh(g::ExtendableGrid; filename::String="")
````

Then the SimplexGrid 'g' is loaded into a gmsh module.
If a string (not "") is passed via 'filename', the mesh is written into this file.

"""
function load_simplexgrid_to_gmsh(g::ExtendableGrid; filename::String="")
    simplexgrid_to_gmshfile(g; filename)
end

"""
````
function load_mixedgrid_to_gmsh(g::ExtendableGrid; filename::String="")
````

Then the ExtendableGrid 'g' is loaded into a gmsh module.
If a string (not "") is passed via 'filename', the mesh is written into this file.

"""
function load_mixedgrid_to_gmsh(g::ExtendableGrid; filename::String="")
    mixedgrid_to_gmshfile(g; filename)
end



#-------------------------------------------------------------------------------------------



"""
````
function gmshfile_to_simplexgrid(filename::String)
````

This function just reads an msh file, and creates a gmsh.model and then calls the 'mod_to_simplexgrid' function
(This function has to be called with an initialized gmsh environment)
This function is called in 'simplexgrid_from_gmsh'
"""
function gmshfile_to_simplexgrid(filename::String; incomplete=false)
    gmsh.open(filename)
    if !incomplete
	    return mod_to_simplexgrid(gmsh.model)
	else
		return incomplete_mod_to_simplexgrid(gmsh.model)
			   
	end
end


"""
````
function gmshfile_to_mixedgrid(filename::String)
````

This function just reads an msh file, and creates a gmsh.model and then calls the 'mod_to_mixedgrid' function
(This function has to be called with an initialized gmsh environment)
This function is called in 'mixedgrid_from_gmsh'
"""
function gmshfile_to_mixedgrid(filename::String)
    gmsh.open(filename)
    return mod_to_mixedgrid(gmsh.model)   
end



"""
````
function simplexgrid_to_gmshfile(grid::ExtendableGrid, filename::String)
````

This function takes a simplexgrid, uses 'grid_to_mod' to create a corresponding gmsh module
Then it writes the module to a file
(this function has to be called with an initialized gmsh environment)
This function is called in 'write_gmsh'
"""
function simplexgrid_to_gmshfile(grid::ExtendableGrid; filename::String="")
	if VERSION>=v"1.9"
        # This possibility is new in 1.9, see
        # https://github.com/JuliaLang/julia/blob/release-1.9/NEWS.md#new-language-features
        gmsh.model = simplexgrid_to_mod(grid)
    else
        simplexgrid_to_mod(grid)
    end
    
    if filename != ""
    	gmsh.write(filename)
    end
end


"""
````
function mixedgrid_to_gmshfile(grid::ExtendableGrid, filename::String)
````

This function takes a mixed grid, uses 'grid_to_mod' to create a corresponding gmsh module
Then it writes the module to a file
(this function has to be called with an initialized gmsh environment)
This function is called in 'write_gmsh'

grid[CellNodes] must be a VariableTargetAdjacency structure
"""
function mixedgrid_to_gmshfile(grid::ExtendableGrid; filename::String="")
	if VERSION>=v"1.9"
        # This possibility is new in 1.9, see
        # https://github.com/JuliaLang/julia/blob/release-1.9/NEWS.md#new-language-features
        gmsh.model = mixedgrid_to_mod(grid)
    else
        mixedgrid_to_mod(grid)
    end
    
    if filename != ""
	    gmsh.write(filename)
    end
end


#---------------------------------------------------------------------------------------------


#=
#old names:!!!!
this file contains the 2 main functions:
a) "mod_to_grid": takes a gmsh.module and creates an ExtendableGrid from it
b) "grid_to_mod": takes an ExtendableGrid and creates a gmsh.module from it

for "mod_to_grid" there also exists "gmshfile_to_grid" which loads the content of an msh file and then calls "mod_to_grid"
for "grid_to_mod" there also exists "grid_to_gmshfile" which loads the content of a grid into a gmsh module and then calls "grid_to_mod"
=#


### general support functions


"""
````
function test_gmsh_init()
````

Very primitive function to test, via a try-catch-block, whether gmsh is already initialized.
If not, it will be initialized.

"""
function test_gmsh_init()
    try
        gmsh.option.setNumber("General.Terminal", 1)
	catch e
		@warn "gmsh may not have been initialized. but is initialized now!"
		gmsh.initialize()
	end
end



"""
````
function take_second(x)
````

x is a list of 2-tuples, with an Int as second entry
an array of the second entries is returned
"""
function take_second(x)
    y = zeros(Int64, length(x))
    for i = 1:length(x)
        _, t = x[i]
        y[i] = t
    end
    return y
end



"""
````
function multiply_indices(indices, n)
````

for n=3:
[i, j, ..., k], 3 -> [3*i-2, 3*i-1, 3*i, 3*j-1, 3*j-2, 3*j, ..., 3*k-2, 3*k-1, 3*k]
in general:
[i, j, ..., k], n -> [n*i-(n-1), n*i-(n-2), ..., n*i, n*j-(n-1), ...]
This function can be used, if you have the indices of cells, and you want to get all their nodes, but the nodes are stored in one list for all cells: 
[node_1_of_cell1, node_2_of_cell1, ... node_n_of_cell1, node_1_of_cell2, ...]
"""
function multiply_indices(indices, n)
    m = length(indices)
    ind_new = zeros(Int64, n * m)
    for i = 1:n
        ind_new[(i-1)*m+1:i*m] = n * indices .- (n - i)
    end
    return sort(ind_new)
end

"""
````
function use_vta(VTA, col_ids, num)
````


If VTA were a matrix, the result would be equivalent to VTA[:, col_ids].
Each column of the VTA contains the nodes of one cell.
"""
function use_vta(VTA, col_ids, num)
	result = zeros(Int64, num*length(col_ids))
	count = 1
	for j in col_ids
		for i=1:num
			result[count] = VTA[i,j]
			count += 1
		end
	end
	return result
end

"""
````
function use_geoms(cellgeoms, ids)
````

If cellgeoms would just be an array/vector, the result would be equivalent to cellgeoms[ids].
"""
function use_geoms(cellgeoms, ids)
	res_cellgeoms = []
	for id in ids
		push!(res_cellgeoms, cellgeoms[id])
	end
	return res_cellgeoms
end




"""
````
function mod_to_grid(model::Module)
````
    
Function that tries to create an (simplex-) ExtendableGrid from a gmsh.model.
Model has to be a gmsh.model.
(This function has to be called with an initialized gmsh environment).
This function is called in 'simplexgrid_from_gmsh'.
    
"""
function mod_to_simplexgrid(model::Module)
    
	dim = model.getDimension()

    node_names, coords, _ = model.mesh.getNodes()
    cell_types, element_names_cells, base_nodes_cells = model.mesh.getElements(dim, -1)
    face_types, element_names_faces, base_nodes_faces = model.mesh.getElements(dim - 1, -1)
    

    #check whether cells are tetrahedrons in 3d or triangles in 2d:
    if dim == 3
        if cell_types[1] != 4
            @warn "3-dim file, but not (just) tetrahedrons as cells!!!"
            return
        end
    elseif dim == 2
        if cell_types[1] != 2
            @warn "2-dim file, but not (just) triangles as cells!!!"
            return
        end
    else
        @warn "dim is neither 3 nor 2"
        return
    end

    #if dim=3, the coordinates (of the nodes) just have to be reshaped
    #for dim=2, the z-coordinate has to be deleted
    if dim == 3
        coords_new = reshape(coords, (3, Int(length(coords) / 3)))
    else
        coords_new = zeros(Float64, 2, Int(length(coords) / 3))
        for i = 1:Int(length(coords) / 3)
            coords_new[:, i] = coords[3*i-2:3*i-1]
        end
    end

    #number of cells
    m = length(element_names_cells[1])
    #the nodes making up the cells is stored in "base_nodes_cells", just in the wrong format
    simplices = reshape(base_nodes_cells[1], (dim + 1, m))

    
    #the physicalnames are currently unused
    cellregion_to_physicalname = Bijection{Int64,String}()
    pgnum_to_physcialname = Dict()
    cr_count = 1

    #the cellregions correspond to the physical groups in which the cells are
    cellregions = ones(Int64, m)
    pgs = take_second(model.getPhysicalGroups(dim))

    for pg in pgs
        name = model.getPhysicalName(dim, pg)
        if length(name) == 0
            name = "$pg"
        end
        pgnum_to_physcialname[pg] = name
        cellregion_to_physicalname[cr_count] = name
        cr_count += 1
    end


    for i = 1:m
        _, _, _, entitytag = model.mesh.getElement(element_names_cells[1][i])
        for pg in pgs
            if entitytag in model.getEntitiesForPhysicalGroup(dim, pg)
                cellregions[i] = cellregion_to_physicalname(pgnum_to_physcialname[pg]) #pg
                break
            end
        end
    end



    # assemble the boundary faces, just reads the faces stored in the msh file
    # for incomplete boundaries, there will be a function to seal them
    
    k = length(element_names_faces[1])
    bfaces = reshape(base_nodes_faces[1], (dim, k))

    # physical groups for bfaces
    bfaceregions = ones(Int64, k)
    pgs = take_second(model.getPhysicalGroups(dim - 1))

    bfaceregion_to_physicalname = Bijection{Int64,String}()
    pgnum_to_physcialname = Dict()
    fr_count = 1

    for pg in pgs
        name = model.getPhysicalName(dim - 1, pg)
        if length(name) == 0
            name = "$pg"
        end
        pgnum_to_physcialname[pg] = name
        bfaceregion_to_physicalname[fr_count] = name
        fr_count += 1
    end

    for i = 1:k
        _, _, _, entitytag = model.mesh.getElement(element_names_faces[1][i])
        for pg in pgs
            if entitytag in model.getEntitiesForPhysicalGroup(dim - 1, pg)
               bfaceregions[i] = bfaceregion_to_physicalname(pgnum_to_physcialname[pg])
               break
            end
        end
    end
    
    
    return simplexgrid(coords_new, simplices, cellregions, bfaces, bfaceregions)
    
end

"""
````
function incomplete_mod_to_simplexgrid(model::Module)
````

Loads an incomplete mesh from a msh file.
Then converts into an ExtendableGrids.
'incomplete' in this context means the boundary is missing.
With the 'ExtendableGrids.seal!(grid::ExtendableGrid)' the boundary can be added.
"""
function incomplete_mod_to_simplexgrid(model::Module)
	dim = gmsh.model.getDimension()

    node_names, coords, _ = gmsh.model.mesh.getNodes()
    cell_types, element_names_cells, base_nodes_cells = gmsh.model.mesh.getElements(dim, -1)
    face_types, element_names_faces, base_nodes_faces = gmsh.model.mesh.getElements(dim - 1, -1)
    

    #check whether cells are tetrahedrons in 3d or triangles in 2d:
    if dim == 3
        if cell_types[1] != 4
            @warn "3-dim file, but not (just) tetrahedrons as cells!!!"
            return
        end
    elseif dim == 2
        if cell_types[1] != 2
            @warn "2-dim file, but not (just) triangles as cells!!!"
            return
        end
    else
        @warn "dim is neither 3 nor 2"
        return
    end

    #if dim=3, the coordinates (of the nodes) just have to be reshaped
    #for dim=2, the z-coordinate has to be deleted
    if dim == 3
        coords_new = reshape(coords, (3, Int(length(coords) / 3)))
    else
        coords_new = zeros(Float64, 2, Int(length(coords) / 3))
        for i = 1:Int(length(coords) / 3)
            coords_new[:, i] = coords[3*i-2:3*i-1]
        end
    end

    #number of cells
    m = length(element_names_cells[1])
    #the nodes making up the cells is stored in "base_nodes_cells", just in the wrong format
    simplices = reshape(base_nodes_cells[1], (dim + 1, m))

	Ti = Int64
	Tc = Float64
	grid = ExtendableGrid{Tc, Ti}()
	grid[Coordinates] = convert(Matrix{Tc}, coords_new)
	grid[CellNodes]   = convert(Matrix{Ti}, simplices)
	grid[CellRegions] = ones(Int64, size(simplices)[2]) #parse(Ti, lines[elems1+1]))

	if dim == 2
		grid[CoordinateSystem] = Cartesian2D
		grid[CellGeometries]   = VectorOfConstants{ElementGeometries,Int}(Triangle2D, num_cells(grid))
	else
		grid[CoordinateSystem] = Cartesian3D
		grid[CellGeometries]   = VectorOfConstants{ElementGeometries,Int}(Tetrahedron3D, num_cells(grid))
	end

		
    
	
	
	return grid
end



"""
````
function mod_to_mixedgrid(model::Module)
````

Function that tries to create a (mixed-) ExtendableGrid from a gmsh.model.
Model has to be a gmsh.model.
(This function has to be called with an initialized gmsh environment).
This function is called in 'mixedgrid_from_gmsh'.
"""
function mod_to_mixedgrid(model::Module)
    elementtypes_nn_2d = Dict(2=>3, 3=>4) #key=elementtype id, val=num nodes
	elementtypes_na_2d = Dict(2=>Triangle2D, 3=>Quadrilateral2D)
	#elementtypes_na_2d = Dict(2=>ExtendableGrids.Triangle2D, 3=>ExtendableGrids.Quadrangle2D)
	#elementtypes2d = [2, 3]
	#numnodes_2d    = [3, 4]


	elementtypes_nn_3d = Dict(4=>4, 5=>8, 6=>6) #, 7=>5) #key=elementtype id, val=num nodes
	elementtypes_na_3d = Dict(4=>Tetrahedron3D, 5=>Hexahedron3D, 6=>Prism3D)
	#elementtypes_na_3d = Dict(4=>ExtendableGrids.Tetrahedron3D, 5=>ExtendableGrids.Hexahedron3D, 6=>ExtendableGrids.Prism3D)
	#elementtypes3d = [4, 5, 6, 7]
	#numnodes_3d    = [4, 8, 6, 5]


    
    dim = model.getDimension()

    node_names, coords, _ = model.mesh.getNodes()
    cell_types, element_names_cells, base_nodes_cells = model.mesh.getElements(dim, -1)
    face_types, element_names_faces, base_nodes_faces = model.mesh.getElements(dim - 1, -1)
    
    Tc = Float64
	Ti = Int64
	grid = ExtendableGrid{Tc, Ti}();
   
	
	VTA = VariableTargetAdjacency()
    cellgeom = 0 #[]
    
    if dim == 3
    	@warn "dim=3 is not supported yet"
    elseif dim == 2
		#@warn "dim=2"
		#cells
		for (ti, cell_type) in enumerate(cell_types)
			nn = elementtypes_nn_2d[cell_type]
			na = elementtypes_na_2d[cell_type]
			
			temp_cellgeom = VectorOfConstants{ElementGeometries,Ti}(na, length(element_names_cells[ti]))
			
			if ti==1
				cellgeom = temp_cellgeom
			else
				cellgeom = vcat(cellgeom, temp_cellgeom)
			end
			
			#@warn nn, na
			for (ci, cell) in enumerate(element_names_cells[ti])
				Base.append!(VTA, base_nodes_cells[ti][nn*(ci-1)+1:nn*ci])
				#push!(cellgeom, na)
			end
		end
		
		
		coords_new = zeros(Float64, 2, Int(length(coords) / 3))
        for i = 1:Int(length(coords) / 3)
            coords_new[:, i] = coords[3*i-2:3*i-1]
        end
        grid[Coordinates] = convert(Matrix{Tc}, coords_new)
		
	end
	
	
	
	
	grid[CellGeometries] = cellgeom 
	grid[CellNodes] = VTA
	grid[CellRegions] = ones(Int64, num_sources(VTA))
	
	grid[BFaceGeometries] = VectorOfConstants{ElementGeometries,Ti}(Edge1D, length(element_names_faces[1]))
	#@warn element_names_faces
	
	k = length(element_names_faces[1])
    grid[BFaceNodes] = convert(Matrix{Ti}, reshape(base_nodes_faces[1], (dim, k)))
	grid[BFaceRegions] = ones(Int64, k)
	
	if dim == 2
		grid[CoordinateSystem] = Cartesian2D
		#grid[CellGeometries]   = VectorOfConstants{ElementGeometries,Int}(Triangle2D, num_cells(grid))
	else
		grid[CoordinateSystem] = Cartesian3D
		#grid[CellGeometries]   = VectorOfConstants{ElementGeometries,Int}(Tetrahedron3D, num_cells(grid))
	end
	
	return grid
end


"""
````
function grid_to_mod(grid::ExtendableGrid)
````

This function writes an ExtendableGrid into a gmsh module.
(This function has to be called with an initialized gmsh environment)
At the moment, this function can only be used from the outside via 'write_gmsh', where the newly created gmsh module is written into a msh file.
"""
function simplexgrid_to_mod(grid::ExtendableGrid)
    # gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    #(fbase,fext)=splitext(filename)
    gmsh.model.add("fbase")


    # formatting the coordinates correctly
    coords = grid[Coordinates]
    dim = size(coords, 1)
    num_nodes = size(coords, 2)
    nodetags = collect(1:num_nodes)
    #types are taken from https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-1-_0028Legacy_0029
    if dim == 3
        coords3d = reshape(coords, 3 * num_nodes) #vec(reshape(coords, (1,3*num_nodes)))
        elementtype_cell = 4 #tetrahedron
        elementtype_face = 2 #triangle
    else
        coords3d = vcat(coords, zeros(Float64, (1, num_nodes)))
        coords3d = reshape(coords3d, 3 * num_nodes)
        elementtype_cell = 2 #triangle
        elementtype_face = 1 #line
    end

    #all nodes are added (to the model and to an entity of dimension 'dim' with tag 1)
    gmsh.model.addDiscreteEntity(dim, 1)
    gmsh.model.mesh.addNodes(dim, 1, nodetags, coords3d, [])
    entitycount = 2

    #Cells
    cells          = grid[CellNodes]
    num_cells      = size(cells, 2)
    cellregions    = grid[CellRegions]
    cm_cellregions = countmap(cellregions)
    celltags0      = collect(num_nodes+1:num_nodes+num_cells) #there is a counter for all elements (nodes, cells, faces) added, since the nodes are already added, we have to start at num_nodes+1
    nodetags0      = reshape(cells, (dim + 1) * num_cells)    #vector of nodes contained in the cells: [node_1_of_cell1, node_2_of_cell1, ... node_n_of_cell1, node_1_of_cell2, ...]

    #we iterate over all cellregions. for each cellregion, we create an entity and add the cells of this cellregion to the entity. Then we add a PhysicalGroup containing the entity
    for cr_dict_entry in cm_cellregions
        gmsh.model.addDiscreteEntity(dim, entitycount)
        cr, num_elements = cr_dict_entry
        array_with_element_ids = findall(z -> z == cr, cellregions)

        #only select those cells which have the right cellregionnumber
        celltags = celltags0[array_with_element_ids]
        nodetags = nodetags0[multiply_indices(array_with_element_ids, dim + 1)]

        gmsh.model.mesh.addElementsByType(entitycount, elementtype_cell, celltags, nodetags)
        gmsh.model.addPhysicalGroup(dim, [entitycount], cr)

        entitycount += 1
    end


    #Faces (basically the same as for the cells)
    bfaces          = grid[BFaceNodes]
    num_bfaces      = size(bfaces, 2)
    bfaceregions    = grid[BFaceRegions]
    cm_bfaceregions = countmap(bfaceregions)
    facetags0       = collect(num_cells+num_nodes+1:num_cells+num_nodes+num_bfaces)
    nodetags0       = reshape(bfaces, dim * num_bfaces)

    for fr_dict_entry in cm_bfaceregions
        gmsh.model.addDiscreteEntity(dim - 1, entitycount)
        fr, num_elements = fr_dict_entry
        array_with_element_ids = findall(z -> z == fr, bfaceregions)

        #only select those cells which have the right bfaceregionnumber
        facetags = facetags0[array_with_element_ids]
        nodetags = nodetags0[multiply_indices(array_with_element_ids, dim)]

        gmsh.model.mesh.addElementsByType(entitycount, elementtype_face, facetags, nodetags)
        gmsh.model.addPhysicalGroup(dim - 1, [entitycount], fr)

        entitycount += 1
    end


    return gmsh.model

end





function mixedgrid_to_mod(grid::ExtendableGrid)

	elementtypes_nn_2d = Dict(2=>3, 3=>4) #key=elementtype id, val=num nodes
	elementtypes_na_2d = Dict(Triangle2D=>2, Quadrilateral2D=>3)
	#elementtypes_na_2d = Dict(2=>ExtendableGrids.Triangle2D, 3=>ExtendableGrids.Quadrangle2D)
	#elementtypes2d = [2, 3]
	#numnodes_2d    = [3, 4]


	elementtypes_nn_3d = Dict(4=>4, 5=>8, 6=>6) #, 7=>5) #key=elementtype id, val=num nodes
	elementtypes_na_3d = Dict(Tetrahedron3D=>4, Hexahedron3D=>5, Prism3D=>6)

    # gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    #(fbase,fext)=splitext(filename)
    gmsh.model.add("fbase")
	
	
	# formatting the coordinates correctly
    coords = grid[Coordinates]
    dim = size(coords, 1)
    num_nodes = size(coords, 2)
    nodetags = collect(1:num_nodes)
    #types are taken from https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-1-_0028Legacy_0029
    if dim == 3
        coords3d = reshape(coords, 3 * num_nodes) #vec(reshape(coords, (1,3*num_nodes)))
        
        elementtypes_nn_face = elementtypes_nn_2d #Dict(2=>3, 3=>4) 
        elementtypes_na_face = elementtypes_na_2d #Dict(Triangle2D=>2, Quadrilateral2D=>3)
		
		elementtypes_nn_cell = elementtypes_nn_3d #Dict(4=>4, 5=>8, 6=>6) 
		elementtypes_na_cell = elementtypes_na_3d #Dict(Tetrahedron3D=>4, Hexahedron3D=>5, Prism3D=>6)
        
        #elementtype_cell = 4 #tetrahedron
        #elementtype_face = 2 #triangle
    else
        coords3d = vcat(coords, zeros(Float64, (1, num_nodes)))
        coords3d = reshape(coords3d, 3 * num_nodes)
        
        elementtypes_nn_face = Dict(1=>2) 
        elementtypes_na_face = Dict(Edge1D=>1)
		
		elementtypes_nn_cell = elementtypes_nn_2d #Dict(2=>3, 3=>4) 
        elementtypes_na_cell = elementtypes_na_2d #Dict(Triangle2D=>2, Quadrilateral2D=>3)
		
        #elementtype_cell = 2 #triangle
        elementtype_face = 1 #line
    end

    #all nodes are added (to the model and to an entity of dimension 'dim' with tag 1)
    gmsh.model.addDiscreteEntity(dim, 1)
    gmsh.model.mesh.addNodes(dim, 1, nodetags, coords3d, [])
    entitycount = 2

    #Cells
    cellgeoms      = grid[CellGeometries]
    VTA            = grid[CellNodes]
    cellregions    = grid[CellRegions]
    cm_cellregions = countmap(cellregions)
    num_cells      = 0
    #@warn cellgeoms
    #@warn cellgeoms[1]
    #@warn cellgeoms[1,2,3]
    
    #we iterate over all cellregions. for each cellregion, we create an entity and add the cells of this cellregion to the entity. Then we add a PhysicalGroup containing the entity
    for cr_dict_entry in cm_cellregions
        gmsh.model.addDiscreteEntity(dim, entitycount)
    	cr, num_elements = cr_dict_entry
        array_with_element_ids = findall(z -> z == cr, cellregions)
		#@warn array_with_element_ids
		types_in_region = use_geoms(cellgeoms, array_with_element_ids) # cellgeoms[array_with_element_ids]
		cm_types = countmap(types_in_region)
		
		for (type_EG_name, num) in cm_types
			type_GMSH_name = elementtypes_na_cell[type_EG_name]
			#@warn type_EG_name, type_GMSH_name
			array_el_ids_in_reg_of_type = findall(z -> z == type_EG_name, use_geoms(cellgeoms, array_with_element_ids))
			
			#@warn "done"
			#@warn array_el_ids_in_reg_of_type
			celltags = array_with_element_ids[array_el_ids_in_reg_of_type] #.+num_nodes
			
			#@warn celltags
			#@warn VTA[celltags]
			
			nodetags = use_vta(VTA, celltags, elementtypes_nn_cell[type_GMSH_name]) #vcat(VTA[celltags]...)
			num_cells += length(celltags)
			celltags .+= num_nodes
			
			gmsh.model.mesh.addElementsByType(entitycount, type_GMSH_name, celltags, nodetags)
		end 

        #gmsh.model.mesh.addElementsByType(entitycount, elementtype_cell, celltags, nodetags)
        gmsh.model.addPhysicalGroup(dim, [entitycount], cr)

        entitycount += 1

    end
    
    #----------------------------------------------------
    #=
    cells          = grid[CellNodes]
    num_cells      = size(cells, 2)
    cellregions    = grid[CellRegions]
    cm_cellregions = countmap(cellregions)
    celltags0      = collect(num_nodes+1:num_nodes+num_cells) #there is a counter for all elements (nodes, cells, faces) added, since the nodes are already added, we have to start at num_nodes+1
    nodetags0      = reshape(cells, (dim + 1) * num_cells)    #vector of nodes contained in the cells: [node_1_of_cell1, node_2_of_cell1, ... node_n_of_cell1, node_1_of_cell2, ...]

    #we iterate over all cellregions. for each cellregion, we create an entity and add the cells of this cellregion to the entity. Then we add a PhysicalGroup containing the entity
    for cr_dict_entry in cm_cellregions
        gmsh.model.addDiscreteEntity(dim, entitycount)
        cr, num_elements = cr_dict_entry
        array_with_element_ids = findall(z -> z == cr, cellregions)

        #only select those cells which have the right cellregionnumber
        celltags = celltags0[array_with_element_ids]
        nodetags = nodetags0[multiply_indices(array_with_element_ids, dim + 1)]

        gmsh.model.mesh.addElementsByType(entitycount, elementtype_cell, celltags, nodetags)
        gmsh.model.addPhysicalGroup(dim, [entitycount], cr)

        entitycount += 1
    end
	=#

    #Faces (basically the same as for the cells)
    bfaces          = grid[BFaceNodes]
    num_bfaces      = size(bfaces, 2)
    bfaceregions    = grid[BFaceRegions]
    cm_bfaceregions = countmap(bfaceregions)
    facetags0       = collect(num_cells+num_nodes+1:num_cells+num_nodes+num_bfaces)
    nodetags0       = reshape(bfaces, dim * num_bfaces)

    for fr_dict_entry in cm_bfaceregions
        gmsh.model.addDiscreteEntity(dim - 1, entitycount)
        fr, num_elements = fr_dict_entry
        array_with_element_ids = findall(z -> z == fr, bfaceregions)

        #only select those cells which have the right bfaceregionnumber
        facetags = facetags0[array_with_element_ids]
        nodetags = nodetags0[multiply_indices(array_with_element_ids, dim)]

        gmsh.model.mesh.addElementsByType(entitycount, elementtype_face, facetags, nodetags)
        gmsh.model.addPhysicalGroup(dim - 1, [entitycount], fr)

        entitycount += 1
    end


    return gmsh.model

end






end
