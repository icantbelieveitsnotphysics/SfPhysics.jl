module SfMatter

using Unitful

import PhysicalConstants.CODATA2018: k_B
import ..SfGeometry: volume

export Material, density, yield_strength, maxwell_boltzmann_peak_speed

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
"""
maxwell_boltzmann_peak_speed(m::Unitful.Mass, T::Unitful.Temperature) = sqrt((2k_B*T)/m) |> u"m/s"

end
