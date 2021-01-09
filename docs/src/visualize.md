# Visualization

Provide visualization  facilities for grids and functions on grids
for various Julia graphics APIs.

In the moment, the following features are available:

- Plotting simplex grids with cell and boundary regions
- Plotting piecewise linear functions on simplex grids

The basic structure is to create a GridVisiualizer `p` which allows to have a 
grid of subplots. The `visualize!` methods then plot into a subplot `p[i,j]`. 
Creation of a GridVisiualizer takes a `Plotter` keyword argument, which allow to
specify plotting module. Supported are:

- From the REPL: GLMakie, PyPlot, VTKView (linux only), Plots (no 3D, no unstructured)
- From Pluto notebooks: GLMakie PyPlot, Plots, WGLMakie (experimental), MeshCat (experimental)

Note   that    functionality   here   mostly   will    be   added   on
necessity.  Eventually, this  will  be factored  out  to a  standalone
package.


## API

```@autodocs
Modules = [ExtendableGrids.GridVisualize]
Pages = ["GridVisualize/gridvisualize.jl"]
```
