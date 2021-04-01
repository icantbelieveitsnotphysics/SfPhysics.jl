module SfGeometry

using Unitful

import Base: length

export Shape, AbstractEllipsoid, Ellipsoid, Spheroid, Sphere, SphericalShell, Cylinder, Cuboid, Cube
export volume, radius, area, length
export unit_x, unit_y, unit_z, sphere_volume, sphere_radius, cylinder_volume, cylinder_radius, cylinder_length, 
	spherical_cap_solid_angle, is_sphere, is_spheroid, is_triaxial, equatorial_radius, polar_radius, 
	cross_sectional_area, spherical_shell_volume, spherical_shell_thickness
	
const unit_x = [1, 0, 0]
const unit_y = [0, 1, 0]
const unit_z = [0, 0, 1]

sphere_volume(r::Unitful.Length) = (4π/3)r^3
sphere_radius(v::Unitful.Volume) = cbrt(3v / 4π)

cylinder_volume(r::Unitful.Length, h::Unitful.Length) = π * r^2 * h
cylinder_radius(v::Unitful.Volume, h::Unitful.Length) = sqrt(v / (π * h))
cylinder_length(v::Unitful.Volume, r::Unitful.Length) = v / (π * r^2)

abstract type Shape end
abstract type AbstractEllipsoid <: Shape end

struct Ellipsoid <: AbstractEllipsoid
	a::Unitful.Length
	b::Unitful.Length
	c::Unitful.Length
end

const Spheroid = Ellipsoid

Spheroid(r_equatorial::Unitful.Length, r_polar::Unitful.Length) = Ellipsoid(r_equatorial, r_equatorial, r_polar)

const Sphere = Ellipsoid

Sphere(radius::Unitful.Length) = Ellipsoid(radius, radius, radius)
Sphere(v::Unitful.Volume) = Sphere(sphere_radius(v))

struct Cylinder <: Shape
	radius::Unitful.Length
	length::Unitful.Length
end

struct Cuboid <: Shape
	length::Unitful.Length
	width::Unitful.Length
	height::Unitful.Length
end

Cube(l::Unitful.Length) = Cuboid(l, l, l)
Cube(v::Unitful.Volume) = Cube(cbrt(v))

struct SphericalShell <: AbstractEllipsoid
	r_inner::Unitful.Length
	r_outer::Unitful.Length
end

SphericalShell(r_inner::Unitful.Length, v::Unitful.Volume) = SphericalShell(r_inner, r_inner + spherical_shell_thickness)

function SphericalShell(s1::AbstractEllipsoid, s2::AbstractEllipsoid)
	if !is_sphere(s1) && is_sphere(s2)
		error("Cannot construct a spherical shell from non-spherical constraints")
	end

	return SphericalShell(min(radius(s1), radius(s2)), max(radius(s1), radius(s2)))
end

volume(x::Ellipsoid) = (4π/3) * x.a * x.b * x.c
volume(x::Cylinder) = cylinder_volume(x.length, x.radius)
volume(x::Cuboid) = x.length * x.width * x.height
volume(x::SphericalShell) = spherical_shell_volume(x.r_inner, x.r_outer)

"""
    radius(x::Ellipsoid)
	
Returns the mean radius of `x`.
"""
radius(x::Ellipsoid) = (x.a + x.b + x.c) / 3

"""
    radius(x::Cylinder)
	
Returns the radius of `x`.
"""
radius(x::Cylinder) = x.radius

"""
    radius(x::SphericalShell)
	
The outer radius of `x`.
"""
radius(x::SphericalShell) = x.r_outer

is_sphere(x::Ellipsoid) = x.a == x.b && x.b == x.c
is_spheroid(x::Ellipsoid) = x.a == x.b || x.a == x.c || x.b == x.c
is_triaxial(x::Ellipsoid) = x.a != x.b && x.a != x.c && x.b != x.c

is_sphere(x::SphericalShell) = true
is_spheroid(x::SphericalShell) = false
is_triaxial(x::SphericalShell) = false

"""
    equatorial_radius(x::Ellipsoid)
	
Returns the equatorial radius of `x`, if `x` is spheroidal.

An exception will be raised if `x` is a triaxial ellipsoid.
"""
function equatorial_radius(x::AbstractEllipsoid)
	if is_sphere(x)
		return radius(x)
	elseif is_triaxial(x)
		throw(ArgumentError("Triaxial ellipsoids do not have a defined equatorial radius"))
	else
		if x.a == x.b || x.a === x.c
			return x.a
		else
			return x.b
		end
	end
end

"""
    polar_radius(x::Ellipsoid)
	
Returns the polar radius of `x`, if `x` is spheroidal.

An exception will be raised if `x` is a triaxial ellipsoid.
"""
function polar_radius(x::AbstractEllipsoid)
	if is_sphere(x)
		return radius(x)
	elseif is_triaxial(x)
		error("Triaxial ellipsoids do not have a defined polar radius")
	else
		if x.a == x.b
			return x.c
		elseif x.a == x.c
			return x.b
		else
			return x.a
		end
	end
end

"""
    area(x::Ellipsoid)
	
Returns the (outer) surface area of `x`, if `x` is spheroidal.

An exception will be raised if `x` is a triaxial ellipsoid.
"""
function area(x::AbstractEllipsoid)
	if is_sphere(x)
		return 4π * radius(x)^2
	elseif is_spheroid(x)
		eq = equatorial_radius(x)
		po = polar_radius(x)
		
		if eq > po
			# oblate
			e = sqrt(1 - po^2 / eq^2)
			return 2π * eq^2 * (1 + (po^2 / (e * eq^2)) * atanh(e))
		else
			# prolate
			e = sqrt(1 - eq^2 / po^2)
			return 2π * eq^2 * (1 + (po / (e * eq)) * asin(e))
		end
	else
		error("Surface area of a triaxial ellipsoid not yet implemented")
	end
end

"""
    area(x::Cuboid)
	
Returns the surface area of `x`.
"""
area(x::Cuboid) = 2 * (x.width * x.length) + 2 * (x.width * x.height) + 2 * (x.height * x.length)

"""
    area(x::Cylinder)
	
Returns the surface area of `x`, including endcaps.
"""
area(x::Cylinder) = 2π * x.radius * x.height + 2π * radius^2

"""
	length(x::Cuboid)
	
Return the length of the longest side of `x`.

Note that this is not necessarily the same as `x.length`, as it is not enforced that `x.length` be `x`'s longest dimension.
"""
length(x::Cuboid) = max(x.length, x.width, x.height)

"""
    length(x::Cylinder)
	
Return the length of `x`.
"""
length(x::Cylinder) = x.length

"""
    length(x::Ellipsoid)

Return the diameter of `x` along its longest axis.
"""
length(x::Ellipsoid) = max(x.a, x.b, x.c) * 2

"""
    length(x::SphericalShell)

Return the diameter of `x` along its longest axis.
"""
length(x::SphericalShell) = radius(x) * 2

"""
    cross_sectional_area(x::Ellipsoid)
	
For spheres and spheroids, return the cross sectional area of `x` in the coronal plane.

Will return an error for triaxial ellipsoids.
"""
function cross_sectional_area(x::AbstractEllipsoid)
	if is_triaxial(x)
		error("Cross section not well defined for triaxial ellipsoids")
	else
		return equatorial_radius(x) * polar_radius(x) * π
	end
end

"""
    spherical_cap_solid_angle(θ)
	
Compute the solid angle of a cone with its apex at the apex of the solid angle and apex angle `2θ`.

See https://en.wikipedia.org/wiki/Solid_angle#Cone,_spherical_cap,_hemisphere
"""
spherical_cap_solid_angle(θ) = 2π * (1 - cos(θ)) * u"sr"

"""
    spherical_cap_solid_angle(h_cone::Unitful.Length, r_cone::Unitful.Length)
	
Compute the solid angle of a cone with its apex at the apex of the solid angle, height `h_cone` and radius `r_cone`.
"""
spherical_cap_solid_angle(h_cone::Unitful.Length, r_cone::Unitful.Length) = spherical_cap_solid_angle(atan(r_cone, h_cone) |> u"m/m")

"""
    spherical_shell_volume(r1::Unitful.Length, r2::Unitful.Length)
	
Return the volume of a spherical shell defined by an inner and outer radius.

Smaller of `r1`, `r2` is used as the inner radius, and the larger is used as the outer radius.
"""
spherical_shell_volume(r1::Unitful.Length, r2::Unitful.Length) = sphere_volume(max(r1, r2)) - sphere_volume(min(r1, r2))

"""
    spherical_shell_thickness(r_inner::Unitful.Length, v::Unitful.Volume)
	
Returns the thickness of a spherical shell with inner radius `r_inner` and volume `v`.
"""
spherical_shell_thickness(r_inner::Unitful.Length, v::Unitful.Volume) = cbrt((3v/4π) + r_inner^3) - r

end
