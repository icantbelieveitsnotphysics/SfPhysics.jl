module SfMatter

using Unitful

import PhysicalConstants.CODATA2018: k_B, g_n
import ..SfGeometry: volume

export Material, density, yield_strength, maxwell_boltzmann_peak_speed, bouyancy

struct Material
	density::Unitful.Density
	yield_strength::Union{Missing, Unitful.Pressure}
end

Material(density::Unitful.Density) = Material(density, missing)

# empty definition for overriding
mass(::Nothing) = nothing

density(x) = mass(x) / volume(x) |> u"kg/m^3"

density(m::Material) = m.density
yield_strength(m::Material) = m.yield_strength

"""
    maxwell_boltzmann_peak_speed(m::Unitful.Mass, T::Unitful.Temperature)
	
Compute the velocity at the peak of the Maxwell-Boltzmann speed distribution for a gas with molecular mass `m` and temperature `T`.

``\\sqrt{\\frac{2k_BT}{m}}``

where ``k_B`` is the Boltzmann constant.
"""
maxwell_boltzmann_peak_speed(m::Unitful.Mass, T::Unitful.Temperature) = sqrt((2k_B*T)/m) |> u"m/s"

"""
	bouyancy(ρ_body::Unitful.Density, ρ_fluid::Unitful.Density; g = g_n)
	
Acceleration due to bouyancy effects experienced by an object with density `ρ_body` in a fluid with density `ρ_fluid`, undergoing acceleration `g` (which defaults to Earth gravity).
"""
bouyancy(ρ_body::Unitful.Density, ρ_fluid::Unitful.Density; g = g_n) = upreferred(g * ρ_fluid / ρ_body)

"""
	bouyancy(ρ_body::Unitful.Density, ρ_fluid::Unitful.Density; g = g_n)
	
Acceleration due to bouyancy effects experienced by an object made of `body` in a fluid composed of `fluid`, undergoing acceleration `g` (which defaults to Earth gravity).
"""
bouyancy(body::Material, fluid::Material; g = g_n) = upreferred(g * density(fluid) / density(body))

"""
    bouyancy(v_body::Unitful.Volume, ρ_fluid::Unitful.Density; g = g_n)
	
Bouyancy force felt by an object of volume `v_body` in a fluid of density `ρ_fluid`, undergoing acceleration `g` (which defaults to Earth gravity).
"""
bouyancy(v_body::Unitful.Volume, ρ_fluid::Unitful.Density; g = g_n) = upreferred(g * ρ_fluid * v_body)

"""
    bouyancy(v_body::Unitful.Volume, ρ_fluid::Unitful.Density; g = g_n)
	
Bouyancy force felt by an object of volume `v_body` in a fluid composed of `fluid`, undergoing acceleration `g` (which defaults to Earth gravity).
"""
bouyancy(v_body::Unitful.Volume, fluid::Material; g = g_n) = upreferred(g * density(fluid) * v_body)

end