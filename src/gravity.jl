module SfGravity

	using PhysicalConstants.CODATA2018: c_0, g_n, G, StefanBoltzmannConstant, ħ, k_B
using Unitful
using UnitfulAstro

export vis_viva, gravity, planetary_mass, planetary_radius, escape_velocity, hill_sphere,
	orbital_velocity, orbital_period, gravitational_binding_energy

function vis_viva(parent_mass::Unitful.Mass, semimajor_axis::Unitful.Length, current_radius::Unitful.Length)
	G*parent_mass*(2/current_radius - 1/semimajor_axis) |> u"m^2/s^2"
end

"""
	gravity(m::Unitful.Mass, r::Unitful.Length)
	
Calculate acceleration due to gravity at distance `r` from point mass `m`

# Example
```julia-repl
julia> gravity(1u"Mearth", 1u"Rearth")
9.798075156340099 m s^-2
```"""
gravity(m::Unitful.Mass, r::Unitful.Length) = G * m / r^2 |> u"m/s/s"
	
"""
	gravity(m1::Unitful.Mass, m2::Unitful.Mass, r::Unitful.Length)
	
Calculate the strength of gravity at a distance `r` between two masses, `m1` and `m2`.

# Example
```julia-repl
julia> gravity(1u"Mearth", 100u"kg", 1u"Rearth")
979.8075156340099 N
```"""
gravity(m1::Unitful.Mass, m2::Unitful.Mass, r::Unitful.Length) = G * m1 * m2 / r^2 |> u"N"
	
"""
	planetary_mass(acc::Unitful.Acceleration, r::Unitful.Length)
	
Calculate the mass of a body required to produce gravitational acceleration `acc` at a distance of `r`.
"""
planetary_mass(acc::Unitful.Acceleration, r::Unitful.Length) = (acc * r^2) / G |> u"kg"
	
"""
	planetary_radius(m::Unitful.Mass, acc::Unitful.Acceleration)
	
Calculate the radius of a planet of mass `m` with surface gravitational acceleration of `acc`.
"""
planetary_radius(m::Unitful.Mass, acc::Unitful.Acceleration) = sqrt((G * m) / acc) |> u"m"

	# kepler's third
	
"""
	orbital_period(m::Unitful.Mass, r::Unitful.Length)
	
Calculate the period of an orbit with semimajor axis `r` about a body with mass `m`.
"""
orbital_period(m::Unitful.Mass, r::Unitful.Length) = sqrt((4π^2 * r^3) / (G * m)) |> u"s"
	
"""
	orbital_radius(m::Unitful.Mass, t::Unitful.Time)
	
Calculate the semimajor axis of an orbit with period `t` about a body with mass `m`.
"""
orbital_radius(m::Unitful.Mass, t::Unitful.Time) = cbrt((t^2 * G * m)/(4π^2)) |> u"m"
	
"""
	orbital_velocity(sma::Unitful.Length, t::Unitful.Time)
	
Calculate the average orbital velocity of a circular orbit of period `t`.

The orbit is assumed to be circular.
"""
orbital_velocity(sma::Unitful.Length, t::Unitful.Time) = 2π * sma / t |> u"km/s"
	
"""
    orbital_velocity(sma::Unitful.Length, t::Unitful.Time, e)
	
Approximate the average orbital velocity of a orbit with eccentricity `e` and period `t`.
"""
function orbital_velocity(sma::Unitful.Length, t::Unitful.Time, e)
	(2π * sma / t) * (1 - e^2 / 4 - 3e^4 / 64 - 5e^6 / 256 - 175e^8 / 16384) |> u"km/s"
end

"""
    orbital_velocity(parent_mass::Unitful.Mass, semimajor_axis::Unitful.Length, current_radius::Unitful.Length)

Approximate the orbital velocity of a body with an eccentric orbit with `semimajor_axis` and current orbital radius `current_radius` about a body of mass `parent_mass`	
"""	
orbital_velocity(parent_mass::Unitful.Mass, semimajor_axis::Unitful.Length, current_radius::Unitful.Length) = sqrt(vis_viva(parent_mass, sem, current_radius)) |> u"km/s"

"""
    orbital_velocity(parent_mass::Unitful.Mass, radius::Unitful.Length)
	
Approximate the orbital velocity of a body with a circular orbit of radius `radius` about a body with mass `parent_mass`.
"""
orbital_velocity(parent_mass::Unitful.Mass, radius::Unitful.Length) = sqrt(G * parent_mass / radius) |> u"km/s"
	
"""
	planetary_mass(orbital_radius::Unitful.Length, orbital_period::Unitful.Time)
	
Calculate the mass of a body orbited by a satellite with semimajor axis `sma` and period `orbital_period`.
"""
planetary_mass(sma::Unitful.Length, orbital_period::Unitful.Time) = (4π^2 * sma^3) / (G * orbital_period^2) |> u"kg"

"""
	escape_velocity(m::Unitful.Mass, r::Unitful.Length)
	
Calculate the velocity required to escape a body with mass `m` from distance `r`.
"""
escape_velocity(m::Unitful.Mass, r::Unitful.Length) = sqrt(2G * m / r) |> u"km/s"
	
"""
	planetary_mass(r::Unitful.Length, v_esc::Unitful.Velocity)
	
Calculate the mass of a body with escape velocity `v_esc` at distance `r`.
"""
planetary_mass(r::Unitful.Length, v_esc::Unitful.Velocity) = v_esc^2 * r / 2G |> u"kg"
	
"""
	planetary_radius(m::Unitful.Mass, v_esc::Unitful.Velocity)
	
Calculate the radius of a body of mass `m` with a surface escape velocity of `v_esc`.
"""
planetary_radius(m::Unitful.Mass, v_esc::Unitful.Velocity) = 2G * m / v_esc^2 |> u"m"

"""
	hill_sphere(m_parent::Unitful.Mass, m::Unitful.Mass, sma::Unitful.Length, e = 0)
	
Calculate the approximate Hill sphere radius at periapse of a body of mass `m` that orbits a body of mass `m_parent` at a distance of `sma` with orbital eccentricity `e`.
"""
function hill_sphere(m_parent::Unitful.Mass, m::Unitful.Mass, sma::Unitful.Length, e = 0)
	sma * (1-e) * cbrt(m / 3m_parent) |> Unitful.unit(sma)
end

"""
	gravitational_binding_energy(m::Unitful.Mass, r::Unitful.Length)
	
Calculate the gravitational binding energy of a body with mass `m` and radius `r`.
"""
gravitational_binding_energy(m::Unitful.Mass, r::Unitful.Length) = (3G*m^2) / 5r |>u"J"

end
