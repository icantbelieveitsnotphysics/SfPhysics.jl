module SfGeometry

using Unitful
import Elliptic: E, F

import Base: length

export Shape, Ellipsoid, TriaxialEllipsoid, Spheroid, Sphere, SphericalShell, Cylinder, Cuboid, Cube
export volume, radius, area, length
export unit_x, unit_y, unit_z, sphere_volume, sphere_radius, cylinder_volume, cylinder_radius, cylinder_length, 
	spherical_cap_solid_angle, is_sphere, is_spheroid, is_triaxial, equatorial_radius, polar_radius, 
	cross_sectional_area, spherical_shell_volume, spherical_shell_thickness, shape, intersecting_circles
	
const unit_x = [1, 0, 0]
const unit_y = [0, 1, 0]
const unit_z = [0, 0, 1]

sphere_volume(r::Unitful.Length) = (4π/3)r^3
sphere_area(r::Unitful.Length) = 4π * r^2
sphere_radius(v::Unitful.Volume) = cbrt(3v / 4π)

cylinder_volume(r::Unitful.Length, h::Unitful.Length) = π * r^2 * h
cylinder_radius(v::Unitful.Volume, h::Unitful.Length) = sqrt(v / (π * h))
cylinder_length(v::Unitful.Volume, r::Unitful.Length) = v / (π * r^2)

abstract type Shape end

shape(x::Shape) = x

abstract type Ellipsoid <: Shape end

struct TriaxialEllipsoid <: Ellipsoid
	a::Unitful.Length
	b::Unitful.Length
	c::Unitful.Length
end

struct Spheroid <: Ellipsoid
	equatorial_radius::Unitful.Length
	polar_radius::Unitful.Length
end

struct Sphere <: Ellipsoid
	radius::Unitful.Length
end

Sphere(v::Unitful.Volume) = Sphere(sphere_radius(v))

Ellipsoid(r::Unitful.Length) = Sphere(r)
Ellipsoid(equatorial_radius::Unitful.Length, polar_radius::Unitful.Length) = 
	equatorial_radius == polar_radius ? Sphere(equatorial_radius) : Spheroid(equatorial_radius, polar_radius)

function Ellipsoid(a::Unitful.Length, b::Unitful.Length, c::Unitful.Length)
	if a == b && b == c
		return Sphere(a)
	elseif a == b || b == c
		return Spheroid(a, c)
	elseif a == c
		return Spheroid(a, b)
	else
		return TriaxialEllipsoid(a, b, c)
	end
end

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

struct SphericalShell <: Ellipsoid
	r_inner::Unitful.Length
	r_outer::Unitful.Length
end

SphericalShell(r_inner::Unitful.Length, v::Unitful.Volume) = SphericalShell(r_inner, r_inner + spherical_shell_thickness)
SphericalShell(s1::Sphere, s2::Sphere) = SphericalShell(min(radius(s1), radius(s2)), max(radius(s1), radius(s2)))

volume(x::TriaxialEllipsoid) = (4π/3) * x.a * x.b * x.c
volume(x::Spheroid) = (4π/3) * x.equatorial_radius^2 * x.polar_radius
volume(x::Sphere) = (4π/3) * x.radius^3
volume(x::SphericalShell) = spherical_shell_volume(x.r_inner, x.r_outer)
volume(x::Cylinder) = cylinder_volume(x.length, x.radius)
volume(x::Cuboid) = x.length * x.width * x.height

"""
    radius(x::TriaxialEllipsoid)
	
Returns the mean radius of `x`.
"""
radius(x::TriaxialEllipsoid) = (x.a + x.b + x.c) / 3

"""
    radius(x::Spheroid)
	
Returns the mean radius of `x`.
"""
radius(x::Spheroid) = (x.equatorial_radius + x.polar_radius) / 2

"""
    radius(x::Spheroid)
	
Radius of `x`.
"""
radius(x::Sphere) = x.radius

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

"""
    equatorial_radius(x::Ellipsoid)
	
Not defined for triaxial TriaxialEllipsoids. Equivalent to the radius for spheres.
"""
equatorial_radius(x::Ellipsoid) = error("Equatorial radius not well defined for general TriaxialEllipsoids.")
equatorial_radius(x::Sphere) = x.radius
equatorial_radius(x::SphericalShell) = x.r_outer
equatorial_radius(x::Spheroid) = x.equatorial_radius

"""
    polar_radius(x::Ellipsoid)
	
Not defined for triaxial TriaxialEllipsoids. Equivalent to the radius for spheres.
"""
polar_radius(x::Ellipsoid) = error("Polar radius not well defined for general TriaxialEllipsoids.")
polar_radius(x::Sphere) = x.radius
polar_radius(x::SphericalShell) = x.r_outer
polar_radius(x::Spheroid) = x.polar_radius

area(x::Sphere) = 4π * x.radius^2
area(x::SphericalShell) = 4π * x.r_outer^2

function area(x::Spheroid)
	eq = equatorial_radius(x)
	po = polar_radius(x)
	
	if eq == po
		return sphere_area(eq)
	elseif eq > po
		# oblate
		e = sqrt(1 - po^2 / eq^2)
		return 2π * eq^2 * (1 + (po^2 / (e * eq^2)) * atanh(e))
	else
		# prolate
		e = sqrt(1 - eq^2 / po^2)
		return 2π * eq^2 * (1 + (po / (e * eq)) * asin(e))
	end
end

function area(x::TriaxialEllipsoid)
	# sort axes ascending
	c, b, a = sort([x.a, x.b, x.c])
	
	if a == b == c
		return sphere_area(a)
	elseif a == b
		return area(Spheroid(a, c))
	elseif a == c
		return area(Spheroid(a, b))
	elseif b == c
		return area(Spheroid(b, a))
	end

	cos_ϕ = c / a
	ϕ = acos(cos_ϕ)
	sin_ϕ = sin(ϕ)
	
	k = sqrt((a^2 * (b^2 - c^2))/(b^2 * (a^2 - c^2)))
	
	e = E(ϕ, k)
	f = F(ϕ, k)
	
	2π * c^2 + ((2π * a * b) / sin_ϕ) * (e * sin_ϕ^2 + f * cos_ϕ^2)
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
area(x::Cylinder) = 2π * x.radius * x.length + 2π * x.radius^2

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
	
The diameter of `x` across its widest part.
"""
length(x::Ellipsoid) = max(equatorial_radius(x), polar_radius(x)) * 2
length(x::TriaxialEllipsoid) = max(x.a, x.b, x.c) * 2

"""
    cross_sectional_area(x::Ellipsoid)
	
The area of the largest cross-section of `x` in the coronal plane.

Not defined for triaxial TriaxialEllipsoids.
"""
cross_sectional_area(x::Ellipsoid) = equatorial_radius(x) * polar_radius(x) * π
cross_sectional_area(::TriaxialEllipsoid) = error("Cross section not well defined for triaxial TriaxialEllipsoids")

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
spherical_cap_solid_angle(h_cone::Unitful.Length, r_cone::Unitful.Length) = spherical_cap_solid_angle(atan(r_cone, h_cone))

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

"""
    intersect_circles(d, r_1, r_2)
    
The area of the intersection of two circles with radii `r_1` ≥ `r_2` and centres separated by `d`

# Example
```jldoctest
julia> intersect_circles(2, 1, 1)
0
julia> intersect_circles(0, 1, 1)
3.141592653589793
"""
function intersect_circles(d, r_1, r_2)
    if r_2 > r_1
        throw(ArgumentError("r_1 must be ≥ r_2"))    
    end
    
    if (d - r_1 - r_2) >= 0
        return 0
    elseif d <= (r_1 - r_2)
        return π * r_1^2
    end
    
    # https://diego.assencio.com/?index=8d6ca3d82151bad815f78addf9b5c1c6
    
    d_1= (r_1^2 - r_2^2 + d^2) / 2d
    d_2 = d - d_1 #(tr^2 - r_1^2 + d^2) / 2d
    
    a_1 = r_1^2 * acos(d_1/r_1) - d_1 * sqrt(r_1^2 - d_1^2)
    a_2 = r_2^2 * acos(d_2/r_2) - d_2 * sqrt(r_2^2 - d_2^2)
    
    a_1 + a_2
end

end
