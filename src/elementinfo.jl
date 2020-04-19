const ElementInfo{T}=Union{Vector{T},VectorOfConstants{T}}

abstract type AbstractElementType end
abstract type Simplex1D <: AbstractElementType end
abstract type Simplex2D <: AbstractElementType end
abstract type Simplex3D <: AbstractElementType end
