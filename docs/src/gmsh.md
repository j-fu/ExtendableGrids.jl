# Gmsh interoperability 
This functionality is in beta stage.
Breaking changes for this API are considered non-breaking for the package.
Therefore, these functions are not exported yet.



## API
These methods become available via a package extension which is loaded together with 
[Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl).
See the general [gmsh documentation](https://gmsh.info/), the [Gmsh reference manual](https://gmsh.info/doc/texinfo/gmsh.html)
and the [Gmsh Julia API source code](https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.jl) for  information.


```@docs
ExtendableGrids.simplexgrid_from_gmsh
ExtendableGrids.simplexgrid_to_gmsh
ExtendableGrids.mixedgrid_from_gmsh
ExtendableGrids.mixedgrid_to_gmsh
```


