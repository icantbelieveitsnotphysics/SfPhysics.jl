module SfThermo

using Unitful, Documenter
	
DocMeta.setdocmeta!(SfThermo, :DocTestSetup, :(using Unitful, ..SfThermo); recursive=true)

import PhysicalConstants.CODATA2018: k_B, N_A, σ
import ..SfUnits: Energy, SpecificEnergy, MolarMass

export maxwell_boltzmann_peak_speed, rms_thermal_velocity, black_body_radiant_flux, black_body_radiated_power,
	ideal_gas_kinetic_energy, ideal_gas_cooling_time, convert_energy, convert_temperature

"""
    maxwell_boltzmann_peak_speed(m::Unitful.Mass, T::Unitful.Temperature)
	
Compute the velocity at the peak of the Maxwell-Boltzmann speed distribution for a gas with molecular mass `m` and temperature `T`.

``\\sqrt{\\frac{2k_BT}{m}}``

where ``k_B`` is the Boltzmann constant.
"""
maxwell_boltzmann_peak_speed(m::Unitful.Mass, T::Unitful.Temperature) = sqrt((2k_B*T)/m) |> u"m/s"

"""
    rms_thermal_velocoty(m::Unitful.Mass, T::Unitful.Temperature)
	
Compute RMS of the Maxwell-Boltzmann speed distribution for a gas with molecular mass `m` and temperature `T`.

``\\sqrt{\\frac{3k_BT}{m}}``

where ``k_B`` is the Boltzmann constant.
"""
rms_thermal_velocity(m::Unitful.Mass, T::Unitful.Temperature) = sqrt((3k_B*T)/m) |> u"m/s"

"""
    black_body_radiant_flux(t::Unitful.Temperature, ϵ::Real = 1, amb::Unitful.Temperature = 0u"K")
	
Radiated power per unit area across all frequencies from a black body with emissivity `ϵ` and temperature `t` 
into an environment with ambient temperature `amb`.

``j^* = \\epsilon \\sigma(T_{hot}^4 - T_{ambient}^4)``

where ``k_B`` is Boltzmann's constant, ``\\sigma`` is the Stefan-Boltzmann constant and ``\\epsilon`` is the radiator emissivity.
"""
black_body_radiant_flux(t::Unitful.Temperature, ϵ::Real = 1, amb::Unitful.Temperature = 0u"K") = upreferred(ϵ * σ * (t^4 - amb^4))

"""
    black_body_radiated_power(t::Unitful.Temperature, a::Unitful.Area, ϵ::Real = 1, amb::Unitful.Temperature = 0u"K")
	
Radiated power across all frequencies from a black body with radiating area `a`, emissivity `ϵ` and temperature `t` 
into an environment with ambient temperature `amb`.

``P = A \\epsilon \\sigma(T_{hot}^4 - T_{ambient}^4)``

where ``k_B`` is Boltzmann's constant, ``\\sigma`` is the Stefan-Boltzmann constant, ``A`` is the radiating area and 
``\\epsilon`` is the radiator emissivity.
"""
black_body_radiated_power(t::Unitful.Temperature, a::Unitful.Area, ϵ::Real = 1, amb::Unitful.Temperature = 0u"K") = upreferred(ϵ * a * σ * (t^4 - amb^4))

"""
    ideal_gas_kinetic_energy(n, t::Unitful.Temperature)
	
Kinetic energy of `n` particles of an ideal gas at temperature `t`.

``E = N \\frac{3}{2}k_BT``

where ``k_B`` is Boltzmann's constant.
"""
ideal_gas_kinetic_energy(n, t::Unitful.Temperature) = 3n * k_B * t / 2 |> u"J"

"""
    ideal_gas_cooling_time(n, t_initial::Unitful.Temperature, t_final::Unitful.Temperature, a::Unitful.Area, ϵ::Real = 1)
	
Time taken for `n` particles of an ideal gas with radiating area `a` and emissivity `ϵ` to cool from a temperature
of `t_initial` to a temperature of `t_final`.

``t_{cooling} = \\frac{Nk_B}{2\\sigma \\epsilon A}\\left[ \\frac{1}{T_{final}^3} - \\frac{1}{T_{initial}^3} \\right]``

where ``k_B`` is Boltzmann's constant, ``\\sigma`` is the Stefan-Boltzmann constant, ``A`` is the radiating area, ``\\epsilon`` is 
the radiator emissivity, and ``T_initial`` and ``T_final`` and the start and end temperatures.
"""
ideal_gas_cooling_time(n, t_initial::Unitful.Temperature, t_final::Unitful.Temperature, a::Unitful.Area, ϵ::Real = 1) = 
	upreferred(((n * k_B) / (2ϵ * σ * a)) * ((1/t_final^3) - (1 / t_initial^3)))

"""
    convert_temperature(e::Energy)
	
Convert a temperature defined as eg. energy per atom `e` to kelvin.

# Example
```jldoctest
julia> convert_temperature(1u"eV")
11604.518121550082 K
```"""
convert_temperature(e::Energy) = e / k_B |> u"K"

"""
    convert_temperature(t::Unitful.Temperature)
	
Convert a temperature `t` in eg. kelvin to eV per atom.

# Example
```jldoctest
julia> convert_temperature(12000u"K")
1.0340799914574215 eV
```"""
convert_temperature(t::Unitful.Temperature) = t * k_B |> u"eV"

"""
    convert_energy(e::SpecificEnergy, mn::MolarMass)
	
Convert specific energy `e` for given material `mn` to eV per atom.

# Example
```jldoctest
julia> convert_energy(1e10u"erg/g", 60u"g/mol")
0.6218561793757305 eV
```"""
convert_energy(e::SpecificEnergy, mn::MolarMass) = e * mn / N_A |> u"eV"

"""
    convert_energy(e::Energy, mn::MolarMass)
	
Convert energy per atom `e` to specific energy for a material with given molar mass `mn`.

# Example
```jldoctest
julia> convert_energy(0.6218561793757305u"eV", 60u"g/mol")
999999.9999999999 J kg^-1
```"""
convert_energy(e::Energy, mn::MolarMass) = (e * N_A) / mn |> u"J/kg"

end
