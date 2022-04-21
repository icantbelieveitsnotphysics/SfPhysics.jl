module SfRelativity

using PhysicalConstants.CODATA2018: c_0, g_n, G
using Unitful
using UnitfulAstro

using ..SfUnits

import ..SfPhysics: kinetic_energy

export lorentz_factor, lorentz_velocity, relativistic_kinetic_energy, relativistic_brachistochrone_transit_time, 
	proper_relativistic_brachistochrone_transit_time, relativistic_velocity, relativistic_delta_v, relativistic_mass_ratio,
    relativistic_brachistochrone_acceleration, rapidity, rapidity_to_velocity

"""
	lorentz_factor(v::Unitful.Velocity)
	
Compute the Lorentz factor associated with the given velocity.

Note that for small values of `v` floating point precision issues can produce substantial
inaccuracies unless velocities are specified as `BigInt` or `BigFloat`.
"""
lorentz_factor(v::Unitful.Velocity) = 1/sqrt(1-(v/c_0)^2) |> u"m/m"
lorentz_factor(t::Unitful.Time, acc::Unitful.Acceleration) = sqrt(1 + (acc*t/c_0)^2) |> u"m/m"
lorentz_factor(d::Unitful.Length, acc::Unitful.Acceleration) = 1 + (acc * d)/c_0^2 |> u"m/m"

"""
    lorentz_velocity(γ)
   
The velocity associated with Lorentz factor `γ`.
"""
lorentz_velocity(γ) = sqrt(1 - (1/γ)^2) * c_0 |> u"c"

rapidity(v::Unitful.Velocity) = atanh(v/c_0)

rapidity_to_velocity(r) = c_0 * tanh(r)

"""
	relativistic_kinetic_energy(m::Unitful.Mass, v::Unitful.Velocity)
	
Compute the kinetic energy of a body with mass `m` travelling at velocity `v`.

Relativistic corrections are included. Note that for small values of `v` floating point precision issues
can produce substantial inaccuracies unless velocities are specified as `BigInt` or `BigFloat`.
"""
relativistic_kinetic_energy(m::SfUnits.Mass, v::Unitful.Velocity) = (m * big(c_0)^2)* (lorentz_factor(v) - 1) |> u"J"

"""
	relativistic_velocity(m::SfUnits.Mass, ke::Unitful.Energy)
	
Velocity of a particle of mass `m` with kinetic energy `ke`, accounting for relativity.
"""
function relativistic_velocity(m::SfUnits.Mass, ke::Unitful.Energy)
	x = 1 + ke / (m * c_0^2)
	return sqrt(c_0^2 - c_0^2 / x^2)
end

"""
    relativistic_transit_time(d::Unitful.Length, a::Unitful.Acceleration)
    
Co-ordinate time (external observer) to travel distance `d` at constant acceleration `a`, accounting for relativistic effects.
"""
relativistic_transit_time(d::Unitful.Length, a::Unitful.Acceleration) = sqrt((d/c_0)^2 + (2 * d / a))

"""
    relativistic_transit_time(d::Unitful.Length, a::Unitful.Acceleration)
    
Proper time to travel distance `d` at constant acceleration `a`, accounting for relativistic effects.
"""
proper_relativistic_transit_time(d::Unitful.Length, a::Unitful.Acceleration) = (c_0 / a) * asinh(a * relativistic_transit_time(d, a) / c)

"""
    relativistic_brachistochrone_transit_time(dist::Unitful.Length, acc::Acceleration)
	
Co-ordinate time (external observer) to travel a distance of `dist` given a constant acceleration of `acc` and a final velocity of 0 relative to the start, accounting for relativistic effects.

``2 \\sqrt { \\left( d \\over 2c \\right)^2 + {d \\over acc} }``
"""
function relativistic_brachistochrone_transit_time(dist::Unitful.Length, acc::Unitful.Acceleration)
	hd = dist/2

	2 * sqrt((hd/c_0)^2 + (2 * hd / acc))
end

"""
    proper_relativistic_brachistochrone_transit_time(d::Unitful.Length, acc::Acceleration)
	
Proper time to travel a distance of `dist` given a constant acceleration of `acc` and a final velocity of 0 relative to the start, accounting for relativistic effects.

``T = \\frac{2c}{g} \\acosh \\left( \\frac{ad}{2c^2} + 1 \\right)``
"""
function proper_relativistic_brachistochrone_transit_time(d::Unitful.Length, acc::Unitful.Acceleration)
	hd = d/2

	2 * (c_0 / acc) * acosh(((acc * hd) / c_0^2) + 1)
end

"""
    relativistic_brachistochrone_acceleration(dist::Unitful.Length, t::Unitful.Time)
    
Required constant acceleration to travel a distance of `dist` in co-ordinate time `t` with a final velocity of 0 relative to the start, accounting for relativistic effects.
"""
function relativistic_brachistochrone_acceleration(d::Unitful.Length, t::Unitful.Time)
	hd = d/2
    ht = t/2
    
    2hd/(ht^2 - (hd/c_0)^2)
end

"""
    relativistic_velocity(t::Unitful.Time, acc::Unitful.Acceleration)
	
Velocity attained by a body undergoing constant acceleration `acc` for co-ordinate time `t`, accounting for relativistic effects.
"""
relativistic_velocity(t::Unitful.Time, acc::Unitful.Acceleration) = (acc*t)/sqrt(1+((acc*t)/c_0)^2|>u"s/s") |>u"m/s"

"""
    relativistic_delta_v(ve::Unitful.Velocity, mr)
	
Delta-V of a rocket with exhaust velocity `ve` and mass ratio `mr` accounting for relativistic effects.
Note that `Δv` cannot exceed `c`... if total delta-V exceeds this value, consider using rapidities to compute it instead.
"""
relativistic_delta_v(ve::Unitful.Velocity, mr) = c_0 * tanh((ve/c_0) * log(mr)) |>u"m/s"

"""
    relativistic_mass_ratio(ve::Unitful.Velocity, Δv::Unitful.Velocity)
    
Mass ratio of a rocket with given delta-V `Δv` and exhaust velocity `ve`, accounting for relativistic effects. 
Note that `Δv` cannot exceed `c`... if total delta-V exceeds this value, consider using rapidities to compute it instead.

``\\exp \\left[ {c \\over v_e} \\tanh^{-1} \\left(\\Delta v \\over c \\right) \\right] = \\frac{m_0}{m_1}``
"""
relativistic_mass_ratio(ve::Unitful.Velocity, Δv::Unitful.Velocity) = exp((c_0 / ve) * atanh(Δv / c_0))
	
"""
    relativistic_mass_ratio(ve::Unitful.Velocity, Δr)
    
Mass ratio of a rocket with given delta-rapidity `Δr` and exhaust velocity `ve`, accounting for relativistic effects. 
As rapidities may be freely added, this is the correct method to use if you're ending up with delta-V values above `c`.

``\\exp \\left[ {c \\over v_e} \\tanh^{-1} \\left(\\Delta v \\over c \\right) \\right] = \\frac{m_0}{m_1}``
"""
relativistic_mass_ratio(ve::Unitful.Velocity, Δr) = exp(Δr * c_0 / ve)

"""
    relativistic_mass_ratio(a::Unitful.Acceleration, T::Unitful.Time)
    
Mass ratio of a photo-drive rocket undergoing acceleration `a` for **proper time** `T`.

``\\exp(\\frac{aT}{c}) = \\frac{m_0}{m_1}``
"""
relativistic_mass_ratio(a::Unitful.Acceleration, T::Unitful.Time) = exp(a*T / c_0) - 1

end
