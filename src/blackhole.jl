module SfBlackHole

using Unitful

import PhysicalConstants.CODATA2018: G, ħ, c_0, k_B

export spaghettification_tensile_force, schwarzschild_radius, hawking_temperature, bekenstein_hawking_luminosity, hawking_evaporation

spaghettification_tensile_force(m_parent::Unitful.Mass, m::Unitful.Mass, l::Unitful.Length, r::Unitful.Length) = (G*m_parent*l*m)/4r^3 |> u"N"

"""
    schwarzschild_radius(m::Unitful.Mass)
	
Radius of a Schwarzschild black hole of mass `m`.
"""
schwarzschild_radius(m::Unitful.Mass) = 2G*m / c_0^2 |> u"m"

"""
    hawking_temperature(m::Unitful.Mass)
	
Black body temperature of a Schwarzschild black hole of mass `m`.
"""
hawking_temperature(m::Unitful.Mass) = (ħ * c_0^3)/(8π * G * m * k_B) |> u"K"

"""
    bekenstein_hawking_luminosity(m::Unitful.Mass)
	
Pure photon luminosity of a Schwarzschild black hole of mass `m`.

This approximation will break down as the temperature of the hole increases to the point where particle emission is possible.
"""
bekenstein_hawking_luminosity(m::Unitful.Mass) = upreferred((ħ * c_0^6)/(15360π * G^2 * m^2))

"""
    hawking_evaporation(m::Unitful.Mass)
	
Time for a Schwarzschild black hole of initial mass `m` to evaporate via Hawking radiation.
"""
hawking_evaporation(m::Unitful.Mass) = upreferred((5120π * G^2 * m^3)/(ħ * c_0^4))

end
