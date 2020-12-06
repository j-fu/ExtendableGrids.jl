# Extendable grid data container for numerical simulations

[![Build status](https://github.com/j-fu/ExtendableGrids.jl/workflows/linux-macos-windows/badge.svg)](https://github.com/j-fu/ExtendableGrids.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://j-fu.github.io/ExtendableGrids.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://j-fu.github.io/ExtendableGrids.jl/dev)


Container structure `ExtendableGrid` with type stable content access and lazy content creation holding data for discretization
grids for finite element and finite volume methods. 
Used by [VoronoiFVM](https://github.com/j-fu/VoronoiFVM.jl) and  [GradientRobustMultiPhysics](https://github.com/chmerdon/GradientRobustMultiPhysics.jl),
a package for novel, gradient robust finite element methods.

Contains:
- Tools to create tensor product grids
- Tools for grid modification
- Visualization methods for various backends
- The differently licensed companion package [SimplexGridFactory](https://github.com/j-fu/SimplexGridFactory.jl) contains an interface which allows to
  create `ExtendableGrid` objects with  [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl) which wraps the Triangle mesh generator
  by J. Shewchuk and [TetGen.jl](https://github.com/JuliaGeometry/TetGen.jl) which wraps the  TetGen mesh generator by H. Si.
  


