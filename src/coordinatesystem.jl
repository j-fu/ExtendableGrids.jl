"""
$(TYPEDEF)
Apex type for coordinate systems
"""
abstract type AbstractCoordinateSystem end

"""
$(TYPEDEF)
1D cartesion coordinate system (unknown `x`)
"""
abstract type Cartesian1D   <: AbstractCoordinateSystem end   

"""
$(TYPEDEF)
2D cartesion coordinate system (unknowns `x,y`)
"""
abstract type Cartesian2D   <: AbstractCoordinateSystem end   

"""
$(TYPEDEF)
2D cartesion coordinate system (unknowns `x,y,z`)
"""
abstract type Cartesian3D   <: AbstractCoordinateSystem end   

"""
$(TYPEDEF)
3D cylindrical coordinate system (unknowns `r,ϕ,z`)
"""
abstract type Cylindrical3D <: AbstractCoordinateSystem end

"""
$(TYPEDEF)
2D cylindrical coordinate system (unknowns `r,z`)
"""
abstract type Cylindrical2D <: AbstractCoordinateSystem end

"""
$(TYPEDEF)
2D polar coordinate system (unknowns `r,ϕ`)
"""
abstract type Polar2D       <: AbstractCoordinateSystem end 

"""
$(TYPEDEF)
1D polar coordinate system (unknown `r`)
"""
abstract type Polar1D       <: AbstractCoordinateSystem end

"""
$(TYPEDEF)
3D spheriacal coordinate system (unknowns `r,ϕ,θ`)
"""
abstract type Spherical3D   <: AbstractCoordinateSystem end

"""
$(TYPEDEF)
1D spheriacal coordinate system (unknown `r`)
"""
abstract type Spherical1D   <: AbstractCoordinateSystem end  #r (integral over ϕ,\theta)

