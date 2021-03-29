module SfObjects

using Unitful, ..SfGeometry, ..SfMatter

import ..SfMatter: mass, density, bouyancy

export mass, bouyancy

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

end
