module SfObjects

using Unitful, ..SfGeometry, ..SfMatter

import ..SfMatter: mass, density, bouyancy
import ..SfGeometry: radius, volume, area, shape

export Object
export mass, bouyancy, area, volume, radius, sphere_of, cube_of

struct Object 
	shape::Shape
	material::Material
end

shape(object::Object) = object.shape

radius(object::Object) = radius(object.shape)
area(object::Object) = area(object.shape)
volume(object::Object) = volume(object.shape)

density(object::Object) = density(object.material)

mass(object::Object) = mass(object.material, object.shape)

"""
    mass(x, s::Shape)
	
Calculate the mass of an object with shape `s` made of `x`, which must have [`density`](@ref)
"""
mass(x, s::Shape) = density(m) * volume(s) |> u"kg"

"""
    bouyancy(s::Shape, fluid; g = g_n)
	
Bouyancy force felt by an object of shape `s` in a medium made of `fluid` (which must have [`density`](@ref)),
undergoing acceleration `g` (which defaults to Earth gravity).
"""
bouyancy(s::Shape, fluid; g = g_n) = bouyancy(volume(s), density(fluid), g)

"""
	sphere_of(material::Material, m::Unitful.Mass)
	
Create a sphere of `material` that has mass `m`.
"""
sphere_of(material::Material, m::Unitful.Mass) = Object(Sphere(upreferred(m / density(material))), material)

"""
	cube_of(material::Material, m::Unitful.Mass)
	
Create a cube of `material` that has mass `m`.
"""
cube_of(material::Material, mass::Unitful.Mass) = Object(Cube(upreferred(mass / density(material))), material)

end
