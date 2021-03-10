module SfGeometry

using Unitful

export Shape, Sphere, Cylinder
export volume, radius
export sphere_volume, sphere_radius, cylinder_volume, cylinder_radius, cylinder_length, spherical_cap_solid_angle

sphere_volume(r::Unitful.Length) = (4π/3)r^3
sphere_radius(v::Unitful.Volume) = cbrt(3v / 4π)

cylinder_volume(r::Unitful.Length, h::Unitful.Length) = π * r^2 * h
cylinder_radius(v::Unitful.Volume, h::Unitful.Length) = sqrt(v / (π * h))
cylinder_length(v::Unitful.Volume, r::Unitful.Length) = v / (π * r^2)

abstract type Shape end

struct Sphere <: Shape
	radius::Unitful.Length
end

struct Cylinder <: Shape
	radius::Unitful.Length
	length::Unitful.Length
end

volume(x::Sphere) = sphere_volume(x.radius)
volume(x::Cylinder) = cylinder_volume(x.length, x.radius)

radius(x::Sphere) = x.radius
radius(x::Cylinder) = x.radius

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
