module SfObjects

using Unitful, ..SfGeometry, ..SfMatter

import ..SfMatter: mass, density, bouyancy
import ..SfGeometry: radius, volume, area

export Object
export mass, bouyancy, area, volume, radius, sphere_of, cube_of

struct Object 
	shape::Shape
	material::Material
end

radius(object::Object) = radius(object.shape)
area(object::Object) = area(object.shape)
volume(object::Object) = volume(object.shape)

density(object::Object) = density(object.material)

mass(object::Object) = mass(object.material, object.shape)

"""
    mass(m::Material, s::Shape)
	
Calculate the mass of an object with shape `s` made of material `m`.
"""
mass(m::Material, s::Shape) = density(m) * volume(s) |> u"kg"

"""
    mass(d::Unitful.Density, s::Shape)
	
Calculate the mass of an object with shape `s` made of a material with density `d`.
"""
mass(d::Unitful.Density, s::Shape) = d * volume(s) |> u"kg"

"""
    bouyancy(s::Shape, ρ_fluid::Unitful.Density; g = g_n)
	
Bouyancy force felt by an object of shape `s` in a fluid of density `ρ_fluid`, undergoing acceleration `g` (which defaults to Earth gravity).
"""
bouyancy(s::Shape, ρ_fluid::Unitful.Density; g = g_n) = bouyancy(volume(s), ρ_fluid, g)

"""
    bouyancy(s::Shape, ρ_fluid::Unitful.Density; g = g_n)
	
Bouyancy force felt by an object of shape `s` in a medium made of `fluid`, undergoing acceleration `g` (which defaults to Earth gravity).
"""
bouyancy(s::Shape, fluid::Material; g = g_n) = bouyancy(volume(s), density(fluid), g)

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
