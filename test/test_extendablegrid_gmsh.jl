using GridVisualize
using GLMakie

path1 = "/.../ExtendableGrids.jl/src/" #enter path of "gsmh_to_extendablegrid.jl" here

path2 = "/.../ExtendableGrids.jl/test/" #enter path of msh files here

include(path1*"gsmh_to_extendablegrid.jl")


#test1: read msh file, create a grid, and visualize the grid
#test2: read msh file, create a grid, write this grid into a file again and then read this file, create a new grid and visualize the new grid
#both functions work for dim=2 and dim=3


function test1(dim)
	if dim == 2
		grid = gmshfile_to_grid(path2*"sto_2d.msh", 2)
	else
		grid = gmshfile_to_grid(path2*"sto_3d.msh", 2)
	end
	
	gridplot(grid; Plotter=GLMakie)
end


function test2(dim)
	if dim == 2
		grid = gmshfile_to_grid(path2*"sto_2d.msh", 2)
	else
		grid = gmshfile_to_grid(path2*"sto_3d.msh", 2)
	end
	
	#or write the grid into a file again:
	grid_to_file(grid, path2*"testwrite.msh")

	#and then read this file again:
	gridcheck = gmshfile_to_grid(path2*"testwrite.msh", 2)

	#and then plot this grid again:
	gridplot(gridcheck; Plotter=GLMakie)
end



