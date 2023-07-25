module ExtendableGridsGmshExt

#if isdefined(Base, :get_extension)
#    ###!!! We do only need gmsh directly, no Gridap stuff
#    import Gmsh: gmsh
#else
#    import ..Gmsh: gmsh
#end

import ExtendableGrids: ExtendableGrid, simplexgrid
import ExtendableGrids: Coordinates, CellNodes, CellRegions, BFaceNodes, BFaceRegions
import ExtendableGrids: simplexgrid_from_gmsh, write_gmsh


#!!! Make a license warning at initialization ? Gmsh is GPL - mention this in the readme.

###??? Do we really need this dependency here ? I would rather like to live without, the more that it seems to make some problems.
using Gmsh: gmsh
using StatsBase: countmap
using Bijections

###!!! I am still thinking about the architecture here. May be we should end up with grid_from_msh

# ExtendableGrids method extensions
function simplexgrid_from_gmsh(filename::String)
	test_gmsh_init()
    gmshfile_to_grid(filename)
end


# ExtendableGrids method extensions
function simplexgrid_from_gmsh(mod::Module)
    test_gmsh_init()
    mod_to_grid(mod)
end

###!!! maybe grid_to_msh, not sure yet
function write_gmsh(filename, g)
    test_gmsh_init()
    grid_to_gmshfile(g, filename)
end


###!!! all the rest is (modulo formatting) the untouched stuff from #29

#=
this file contains the 2 main functions:
a) "mod_to_grid": takes a gmsh.module and creates an ExtendableGrid from it
b) "grid_to_mod": takes an ExtendableGrid and creates a gmsh.module from it

for "mod_to_grid" there also exists "gmshfile_to_grid" which loads the content of an msh file and then calls "mod_to_grid"
for "grid_to_mod" there also exists "grid_to_gmshfile" which loads the content of a grid into a gmsh module and then calls "grid_to_mod"
=#


### general support functions


function test_gmsh_init()

    try
    	#@warn "try to use gsmh"
	    #gmsh.initialize()
	    gmsh.option.setNumber("General.Terminal", 1)
	catch e
		@warn "gmsh may not have been initialized. but is initialized now!"
		gmsh.initialize()
	end
	#path = "/home/johannes/Nextcloud/Documents/Uni/VIII/WIAS/ExtendableGrids.jl_in_src/test/"
	    
    #grid1 = ExtendableGrids.simplexgrid_from_gmsh(path*"sto_2d.msh")
    
    #ExtendableGrids.grid_to_gmshfile(grid1, path*"testfile.msh")
    
    #grid2 = ExtendableGrids.simplexgrid_from_gmsh(path*"testfile.msh")
    
    #gmsh.finalize()
    
    #@test seemingly_equal2(grid1, grid2;confidence=:low)
    #@test seemingly_equal2(grid1, grid2;confidence=:medium)
end



"""
    take_second(x)
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
    multiply_indices(indices, n)
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
    mod_to_grid(model)
    function that tries to create an (simplex-) extendable grid from an gmsh.model
    model has to be a gmsh.model
    (this function has to be called with an initialized gmsh environment)
    
"""
function mod_to_grid(model::Module)
    
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
    gmshfile_to_grid(filename::String)
    this function just reads an msh file, and creates a gmsh.model and then calls the 'mod_to_grid' function
    (this function has to be called with an initialized gmsh environment)
"""
function gmshfile_to_grid(filename::String)
   
    gmsh.open(filename)

    grid = mod_to_grid(gmsh.model)

    return grid
end


"""
    grid_to_mod(grid::ExtendableGrid)
    this function writes an ExtendableGrid into a gmsh module
    (this function has to be called with an initialized gmsh environment)
"""
function grid_to_mod(grid::ExtendableGrid)
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


"""
    grid_to_gmshfile(grid::ExtendableGrid, filename::String)
    this function takes a grid, uses 'grid_to_mod' to create a corresponding gmsh module
    then it writes the module to a file
    (this function has to be called with an initialized gmsh environment)
"""
function grid_to_gmshfile(grid::ExtendableGrid, filename::String)
    
    gmsh.model = grid_to_mod(grid)

    gmsh.write(filename)
    
end



end
