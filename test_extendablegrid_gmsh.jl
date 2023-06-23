using GridVisualize
using GLMakie

path = "/..." #enter path here

include(path*"gsmh_to_extendablegrid.jl")

#load the grid from the file
grid3d = gmshfile_to_grid(path*"sto3d_pg_0.1.msh", 2)

testversion = false 
#if testversion=true:  the grid read from the file is visualized
#if testversion=false: the grid is written to a file again and then this file is read and plot again

if testversion
	#plot the grid
	gridplot(grid3d; Plotter=GLMakie)

else
	#or write the grid into a file again:
	grid_to_file(grid3d)

	#and then read this file again:
	gridcheck = gmshfile_to_grid("testwrite.msh", 2)

	#and then plot this grid again:
	gridplot(gridcheck; Plotter=GLMakie)
end


