# Extendable grid data container for numerical simulations

[![Build Status](https://img.shields.io/travis/j-fu/ExtendableGrids.jl/master.svg?label=Linux+MacOSX+Windows)](https://travis-ci.com/j-fu/ExtendableGrids.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://j-fu.github.io/ExtendableGrids.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://j-fu.github.io/ExtendableGrids.jl/dev)


Container structure with type stable content access and lazy content creation holding data for discretization
grids for finite element and finite volume methods. 
Used by [VoronoiFVM](https://github.com/j-fu/VoronoiFVM.jl) and some simulation codes under development.

Contains:
- interface to [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl)
- Tools to create tensor product grids
- Tools for grid modification
- Visualization methods for various backends


