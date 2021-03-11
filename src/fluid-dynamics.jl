module SfFluidDynamics

using Unitful

import PhysicalConstants.CODATA2018: g_n

import ..SfPhysics: Angle, to_angle
import ..SfPhysics: projectile_displacement, projectile_velocity, projectile_flight_time, 
	projectile_peak_displacement, projectile_range, projectile_angle

export dynamic_viscosity, drag_force, reynolds_number

const Viscosity{T} = Unitful.AbstractQuantity{T, Unitful.𝐌 * Unitful.𝐋^-1 *  Unitful.𝐓^-1, typeof(u"Pa*s")}

dynamic_viscosity(drag_coefficient::Real, ρ_fluid::Unitful.Density, projectile_area::Unitful.Area, projectile_mass::Unitful.Mass) =
	(drag_coefficient * ρ_fluid * projectile_area) / 2projectile_mass
	
drag_force(drag_coefficient::Real, ρ_fluid::Unitful.Density, reference_area::Unitful.Area, v::Unitful.Velocity) = 
	0.5ρ_fluid * v^2 * drag_coefficient * reference_area

reynolds_number(v::Unitful.Velocity, ρ_fluid::Unitful.Density, l::Unitful.Length, μ::Viscosity) = (v * ρ_fluid * l) / μ |>u"s/s"
	
end
