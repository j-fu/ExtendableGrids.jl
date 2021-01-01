# Plotting

Provide plotting facilities for various Julia graphics APIs.

In the moment, the following features are available:

- Plotting simplex grids with cell and boundary regions
- Plotting piecewise linear functions on simplex grids

The basic structure is to create a PlotContext `p` which allows to have a 
grid of subplots. The `plot!` methods then plot into a subplot `p[i,j]`. 
Creation of a plot context takes a `Plotter` keyword argument, which allow to
specify plotting module. Supported are:

- From the REPL: GLMakie, PyPlot, VTKView (linux only)
- From Pluto notebooks: GLMakie PyPlot, WGLMakie (experimental), MeshCat (experimental)

Note   that    functionality   here   mostly   will    be   added   on
necessity.  Eventually, this  will  be factored  out  to a  standalone
package.


## API

```@autodocs
Modules = [ExtendableGrids]
Pages = ["plot.jl"]
```
