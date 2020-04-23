
abstract type AbstractCoordinateSystem end
abstract type Cartesian1D   <: AbstractCoordinateSystem end    #x 
abstract type Cartesian2D   <: AbstractCoordinateSystem end    #x,y
abstract type Cartesian3D   <: AbstractCoordinateSystem end    #x,y,z
abstract type Cylindrical2D <: AbstractCoordinateSystem end  #r,z
abstract type Cylindrical3D <: AbstractCoordinateSystem end  #r,ϕ,z
abstract type Polar2D       <: AbstractCoordinateSystem end  #r,ϕ
abstract type Polar1D       <: AbstractCoordinateSystem end  #r (integral over ϕ)
abstract type Spherical3D   <: AbstractCoordinateSystem end  #r,ϕ,θ
abstract type Spherical1D   <: AbstractCoordinateSystem end  #r (integral over ϕ,\theta)

