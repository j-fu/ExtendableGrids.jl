module ExtendableGridsGmshExt

if isdefined(Base, :get_extension)
    ###!!! We do only need gmsh directly, no Gridap stuff
    import Gmsh: gmsh
else
    import ..Gmsh: gmsh
end

import ExtendableGrids: ExtendableGrid, simplexgrid
import ExtendableGrids: Coordinates, CellNodes, CellRegions, BFaceNodes, BFaceRegions
import ExtendableGrids: simplexgrid_from_gmsh, write_gmsh

#!!! Make a license warning at initialization ? Gmsh is GPL - mention this in the readme.

###??? Do we really need this dependency here ? I would rather like to live without, the more that it seems to make some problems.
using StatsBase: countmap
using Bijections

###!!! I am still thinking about the architecture here. May be we should end up with grid_from_msh

# ExtendableGrids method extensions
function simplexgrid_from_gmsh(filename::String)
    gmshfile_to_grid(filename, 0)
end


# ExtendableGrids method extensions
function simplexgrid_from_gmsh(mod::Module)
    mod_to_grid(mod, 0)
end

###!!! maybe grid_to_msh, not sure yet
function write_gmsh(filename, g)
    grid_to_file(g, filename)
end

###!!! all the rest is (modulo formatting) the untouched stuff from #29
#=
this file contains the 2 main functions:
a) "mod_to_grid": takes a gmsh.module and creates an ExtendableGrid from it
b) "grid_to_file": takes an ExtendableGrid and writes the grid into an msh file

for "mod_to_grid" there also exists "gmshfile_to_grid" which loads the content of an msh file and then calls "mod_to_grid"
=#


### general support functions



function take_second(x)
    y = zeros(Int64, length(x))
    for i = 1:length(x)
        _, t = x[i]
        y[i] = t
    end
    return y
end

function set_manual(manualset, face_types, dim)
    if manualset == 2
        if dim == 3
            if length(face_types) == 0 || face_types[1] != 2
                @warn "3-dim file, but not (just) triangles as faces!!!"
                @warn "trying to collect them manually"
                return true
            else
                @warn "manual = false"
                return false
            end

        else #dim == 2
            if length(face_types) == 0 || face_types[1] != 1
                @warn "2-dim file, but not (just) lines as faces!!!"
                @warn "trying to collect them manually"
                return true
            else
                @warn "manual = false"
                return false
            end
        end
    elseif manualset == 0
        @warn "setting: manual = false"
        return false
    elseif manualset == 1
        @warn "setting: manual = true"
        return true
    end
end

function cut_it(simplices, coords, xlim)
    simps = []
    for i = 1:size(simplices, 2)
        is = simplices[:, i]
        #@warn is
        if sum(coords[1, is]) / 4 > xlim
            push!(simps, is)
        end
    end
    return hcat(simps...)
end

function multiply_indices(indices, n)
    m = length(indices)
    ind_new = zeros(Int64, n * m)
    for i = 1:n
        ind_new[(i-1)*m+1:i*m] = n * indices .- (n - i)
    end
    return sort(ind_new)
end


### if the file does not contain all the boundary faces for the elements, then they are assembled manually
###??? not sure if we really will need this. In the moment, the tool stack assumes boundaries.
###??? I tend to hand over the responsibility to the user, and in the extreme case allow grids without boundaries
###??? eventually there may be a call "make_boundary!" or "seal!" which does this.

###!!! Please ensure that orientationwise etc this is compatible to the rest.
function faces_of_ndim_simplex(x, dim, nn)
    # nn = number of total nodes
    # for a given array x with n entries. for n=4: x1, x2, x3, x4
    # i want all subsets with length n-1 (and distinct entries)
    if dim == 3 # =faces_of_tetrahedron
        y = [
            encode(sort([x[1], x[2], x[3]]), nn),
            encode(sort([x[1], x[2], x[4]]), nn),
            encode(sort([x[1], x[3], x[4]]), nn),
            encode(sort([x[2], x[3], x[4]]), nn),
        ]
    elseif dim == 2
        y = [
            encode(sort([x[1], x[2]]), nn),
            encode(sort([x[1], x[3]]), nn),
            encode(sort([x[2], x[3]]), nn),
        ]
    end
    return y

end

function encode(x, nn)
    y = 0
    for i = 1:length(x)
        y += x[i] * nn^(i - 1)
    end
    return y #x[1]+nn*x[2]+nn*nn*x[3]
end

function decode(y, nn, dim)
    x = zeros(Int64, dim)
    x[1] = y % nn
    x[2] = Int(round(y / nn - 0.5) % nn)
    if dim == 3
        x[3] = Int(round(y / nn^2 - 0.5))
    end
    return x #[y%nn, Int(round(y/nn-0.5)%nn), Int(round(y/nn^2-0.5))]
end

#!!!
#!!! I think this could really part of the main extendablegrids module.
#!!! 
function assemble_bfaces(simplices, dim, nn)
    m = size(simplices, 2)
    poss_faces = zeros(Int64, (dim + 1) * m)
    for i = 1:m
        poss_faces[(dim+1)*i-dim:(dim+1)*i] =
            faces_of_ndim_simplex(simplices[:, i], dim, nn)
    end
    dict = countmap(poss_faces)

    k = 0
    for d in dict
        (a, b) = d
        if b == 1
            k += 1
            #push!(unicats, a)
        end
    end

    bfaces = zeros(Int64, (dim, k))
    k = 1
    for d in dict
        (a, b) = d
        if b == 1
            bfaces[:, k] = decode(a, nn, dim)
            k += 1
        end
    end

    return bfaces
end


### function that tries to create an extendable grid from an gmsh.model
###!!! please add a real docstring here.
###!!! Ah manualset is the "seal!" call. I think this should be an extra call outside of this ext. 
function mod_to_grid(model::Module, manualset::Int64)
    #model: gmsh.model
    #(this function has to be called, before the gmsh environment is finalized)
    #manual: 
    #        false->faces are taken from "getElements(dim-1)"
    #        true ->faces are built manually from the cells (expensive)
    #manualset: 0 -> manual = false 
    # 			1 -> manual = true
    # 			2 -> if there are any elements with dim-1, then manual=false,
    # 				 else: manual=true

    dim = model.getDimension()

    node_names, coords, _ = model.mesh.getNodes()
    A = model.mesh.getElements(dim, -1)
    cell_types, element_names_cells, base_nodes_cells = A
    B = model.mesh.getElements(dim - 1, -1)
    face_types, element_names_faces, base_nodes_faces = B

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

    #decide how to set manual, based on manualset
    manual = set_manual(manualset, face_types, dim)
    m = length(element_names_cells[1])
    #the nodes making up the cells is stored in "base_nodes_cells", just in the wrong format
    simplices = reshape(base_nodes_cells[1], (dim + 1, m))

    #look at a small part of the thing
    #simplices   =  cut_it(simplices, coords_new) #simplices[:,1:1000]
    #cellregions = ones(Int64, size(simplices, 2))

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



    #assemble the boundary faces)
    if manual
        bfaces = assemble_bfaces(simplices, dim, Int(length(coords) / 3))
        bfaceregions = ones(Int64, size(bfaces, 2))
    else
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
    end



    return simplexgrid(coords_new, simplices, cellregions, bfaces, bfaceregions)
end

### this function just reads an msh file, and creates a gmsh.model and then calls the 'mod_to_grid' function

function gmshfile_to_grid(filename::String, manualset::Int64)
    #filename of msh file
    #manual: 
    #        false->faces are taken from "getElements(dim-1)"
    #        true ->faces are built manually from the cells (expensive)
    #manualset: 0 -> manual = false 
    # 			1 -> manual = true
    # 			2 -> if there are any elements with dim-1, then manual=false,
    # 				 else: manual=true
    # gmsh.initialize() #!!! this should be somewehere else ?
    gmsh.open(filename)

    grid = mod_to_grid(gmsh.model, manualset)

    #gmsh.finalize() #!!! this should be somewehere else

    return grid
end


### this function writes an ExtendableGrid into a gmsh file with the name "filename"
### !!! split this - make an extra ExtendableGridToGmsh call
function grid_to_file(grid::ExtendableGrid, filename::String)
    # gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    (fbase,fext)=splitext(filename)
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
    cells = grid[CellNodes]
    num_cells = size(cells, 2)
    cellregions = grid[CellRegions]
    cm_cellregions = countmap(cellregions)
    celltags0 = collect(num_nodes+1:num_nodes+num_cells)
    nodetags0 = reshape(cells, (dim + 1) * num_cells)

    for cr_dict_entry in cm_cellregions
        #we iterate over all cellregions. for each cellregion, we create an entity and add the cells of this cellregion to the entity. Then we add a PhysicalGroup containing the entity
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


    #Faces
    bfaces = grid[BFaceNodes]
    num_bfaces = size(bfaces, 2)
    bfaceregions = grid[BFaceRegions]
    cm_bfaceregions = countmap(bfaceregions)
    facetags0 = collect(num_cells+num_nodes+1:num_cells+num_nodes+num_bfaces)
    nodetags0 = reshape(bfaces, dim * num_bfaces)

    #@warn "cm bfaces : $cm_bfaceregions"
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


    gmsh.write(filename)

    # gmsh.finalize()

end




end

