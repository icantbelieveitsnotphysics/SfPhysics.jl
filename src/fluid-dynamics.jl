module SfFluidDynamics

using Unitful

import PhysicalConstants.CODATA2018: g_n

import ..SfPhysics: Angle, to_angle
import ..SfPhysics: projectile_displacement, projectile_velocity, projectile_flight_time, 
	projectile_peak_displacement, projectile_range, projectile_angle

export drag_force, reynolds_number

const Viscosity{T} = Unitful.AbstractQuantity{T, Unitful.𝐌 * Unitful.𝐋^-1 *  Unitful.𝐓^-1, typeof(u"Pa*s")}
	
"""
    drag_force(drag_coefficient::Real, ρ_fluid::Unitful.Density, reference_area::Unitful.Area, v::Unitful.Velocity)
	
Drag acting on a body with `drag_coefficient` and `reference_area` travelling at `v` through a fluid with density `ρ_fluid`.
"""
drag_force(drag_coefficient::Real, ρ_fluid::Unitful.Density, reference_area::Unitful.Area, v::Unitful.Velocity) = 
	upreferred(0.5ρ_fluid * v^2 * drag_coefficient * reference_area)

"""
    reynolds_number(v::Unitful.Velocity, ρ_fluid::Unitful.Density, l::Unitful.Length, μ::Viscosity)
	
Reynolds number of a body with characteristic dimension `l` travelling at velocity `v` through a fluid of density `ρ_fluid` and viscosity `μ`.
"""
reynolds_number(v::Unitful.Velocity, ρ_fluid::Unitful.Density, l::Unitful.Length, μ::Viscosity) = (v * ρ_fluid * l) / μ |>u"s/s"
	
end
