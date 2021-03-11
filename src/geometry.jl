module SfGeometry

using Unitful

export Shape, Ellipsoid, Spheroid, Sphere, Cylinder, Cuboid, Cube
export volume, radius, area
export sphere_volume, sphere_radius, cylinder_volume, cylinder_radius, cylinder_length, spherical_cap_solid_angle

sphere_volume(r::Unitful.Length) = (4π/3)r^3
sphere_radius(v::Unitful.Volume) = cbrt(3v / 4π)

cylinder_volume(r::Unitful.Length, h::Unitful.Length) = π * r^2 * h
cylinder_radius(v::Unitful.Volume, h::Unitful.Length) = sqrt(v / (π * h))
cylinder_length(v::Unitful.Volume, r::Unitful.Length) = v / (π * r^2)

abstract type Shape end

struct Ellipsoid <: Shape
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

volume(x::Ellipsoid) = (4π/3) * x.a * x.b * x.c
volume(x::Cylinder) = cylinder_volume(x.length, x.radius)
volume(x::Cuboid) = x.length * x.width * x.height

radius(x::Ellipsoid) = (x.a + x.b + x.c) / 3
radius(x::Cylinder) = x.radius

is_sphere(x::Ellipsoid) = x.a == x.b && x.b == x.c
is_spheroid(x::Ellipsoid) = x.a == x.b || x.a == x.c || x.b == x.c
is_triaxial(x::Ellipsoid) = x.a != x.b && x.a != x.c && x.b != x.c

function equatorial_radius(x::Ellipsoid)
	if is_sphere(x)
		return x.a
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

function polar_radius(x::Ellipsoid)
	if is_sphere(x)
		return x.a
	elseif is_triaxial(x)
		throw(ArgumentError("Triaxial ellipsoids do not have a defined polar radius"))
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

function area(x::Ellipsoid)
	if is_sphere(x)
		return 4π * x.a^2
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
		throw(ArgumentError("Surface area of a triaxial ellipsoid not yet implemented"))
	end
end

area(x::Cuboid) = 2 * (x.width * x.length) + 2 * (x.width * x.height) + 2 * (x.height * x.length)

area(x::Cylinder) = 2π * x.radius * x.height + 2π * radius^2

"""
    spherical_cap_solid_angle(θ)
	
Compute the solid angle of a cone with its apex at the apex of the solid angle and apex angle `2θ`.

See https://en.wikipedia.org/wiki/Solid_angle#Cone,_spherical_cap,_hemisphere
"""
spherical_cap_solid_angle(θ) = 2π * (1 - cos(θ))

"""
    spherical_cap_solid_angle(h_cone::Unitful.Length, r_cone::Unitful.Length)
	
Compute the solid angle of a cone with its apex at the apex of the solid angle, height `h_cone` and radius `r_cone`.
"""
spherical_cap_solid_angle(h_cone::Unitful.Length, r_cone::Unitful.Length) = spherical_cap_solid_angle(atan(r_cone, h_cone) |> u"m/m")

end
