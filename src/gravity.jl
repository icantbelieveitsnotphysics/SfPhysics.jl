module SfGravity

using PhysicalConstants.CODATA2018: c_0, g_n, G, StefanBoltzmannConstant, ħ, k_B
using Unitful, UnitfulAstro, Documenter

export vis_viva, gravity, planetary_mass, planetary_radius, escape_velocity, hill_sphere,
	orbital_velocity, orbital_period, gravitational_binding_energy, roche_limit, gravity_tug_mass,
	barycentric_distance, tidal_acceleration,
	radial_orbit_time, radial_orbit_displacement
	
DocMeta.setdocmeta!(SfGravity, :DocTestSetup, :(using Unitful, UnitfulAstro, ..SfGravity); recursive=true)

"""
    function vis_viva(parent_mass::Unitful.Mass, semimajor_axis::Unitful.Length, current_radius::Unitful.Length)
	
Return the square of the relative velocities of a body with negligible mass and orbit with `semimajor_axis` and `current_radius` from a central body of `parent_mass`.

``v^2 = GM\\left( \\frac{2}{r} - \\frac{1}{a} \\right)``

where ``G`` is the gravitational constant, ``M`` is the mass of the central body, ``r`` is the current separation of the bodies and ``a`` is the semi-major axis of the orbit.
"""
function vis_viva(parent_mass::Unitful.Mass, semimajor_axis::Unitful.Length, current_radius::Unitful.Length)
	G*parent_mass*(2/current_radius - 1/semimajor_axis) |> u"m^2/s^2"
end

"""
	gravity(m::Unitful.Mass, r::Unitful.Length)
	
Calculate acceleration due to gravity at distance `r` from point mass `m`

# Example
```jldoctest
julia> gravity(1u"Mearth", 1u"Rearth")
9.798398133669465 m s^-2
```"""
gravity(m::Unitful.Mass, r::Unitful.Length) = G * m / r^2 |> u"m/s/s"
	
"""
	gravity(m1::Unitful.Mass, m2::Unitful.Mass, r::Unitful.Length)
	
Calculate the strength of gravity at a distance `r` between two masses, `m1` and `m2`.

# Example
```jldoctest
julia> gravity(1u"Mearth", 100u"kg", 1u"Rearth")
979.8398133669466 N
```"""
gravity(m1::Unitful.Mass, m2::Unitful.Mass, r::Unitful.Length) = G * m1 * m2 / r^2 |> u"N"
	
"""
	planetary_mass(g::Unitful.Acceleration, r::Unitful.Length)
	
Calculate the mass of a body required to produce gravitational acceleration `g` at a distance of `r`.
"""
planetary_mass(g::Unitful.Acceleration, r::Unitful.Length) = (g * r^2) / G |> u"kg"

"""
	planetary_mass(ve::Unitful.Velocity, r::Unitful.Length)
	
Calculate the mass of a body required to produce escape velocity `ve` at a distance of `r`.
"""
planetary_mass(ve::Unitful.Velocity, r::Unitful.Length) = r * ve^2 / 2G |> u"kg"
	
"""
	planetary_radius(m::Unitful.Mass, g::Unitful.Acceleration)
	
Radius of a spherical planet of mass `m` with surface gravitational acceleration of `g`.
"""
planetary_radius(m::Unitful.Mass, g::Unitful.Acceleration) = sqrt((G * m) / g) |> u"m"

"""
	planetary_radius(m::Unitful.Mass, g::Unitful.Acceleration)
	
Radius of a spherical planet of mass `m` with surface escape velocity `ve`.
"""
planetary_radius(m::Unitful.Mass, ve::Unitful.Velocity) = 2G * m / ve^2 |> u"m"

"""
	planetary_radius(ρ::Unitful.Density, g::Unitful.Acceleration)
	
Radius of a spherical planet of density `ρ` with surface gravitational acceleration of `g`.
"""
planetary_radius(ρ::Unitful.Density, g::Unitful.Acceleration) = 3g / (4G * π * ρ) |> u"m"

"""
	planetary_radius(ρ::Unitful.Density, ve::Unitful.Velocity)
	
Radius of a spherical planet of density `ρ` with surface escape velocity `ve`.
"""
planetary_radius(ρ::Unitful.Density, ve::Unitful.Velocity) = sqrt(3ve^2 / (8G * π * ρ)) |> u"m"

# kepler's third
	
"""
	orbital_period(m::Unitful.Mass, r::Unitful.Length)
	
Calculate the period of an orbit with semimajor axis `r` about a body with mass `m`.

``P = \\sqrt{\\frac{4π^2r^3}{Gm}}``
"""
orbital_period(m::Unitful.Mass, r::Unitful.Length) = sqrt((4π^2 * r^3) / (G * m)) |> u"s"
	
"""
	orbital_radius(m::Unitful.Mass, t::Unitful.Time)
	
Calculate the semimajor axis of an orbit with period `t` about a body with mass `m`.

``a = \\sqrt[3]{\\frac{t^2Gm}{4π^2}}``
"""
orbital_radius(m::Unitful.Mass, t::Unitful.Time) = cbrt((t^2 * G * m)/(4π^2)) |> u"m"

"""
    orbital_radius(m::Unitful.Mass, v::Unitful.Velocity)
	
Approximate the radius of a circular orbit with velocity `v` about a body with mass `m`.
"""
orbital_radius(m::Unitful.Mass, v::Unitful.Velocity) = (G * m) / v^2 |> u"m"
	
"""
	orbital_velocity(sma::Unitful.Length, t::Unitful.Time)
	
Calculate the average orbital velocity of a circular orbit of period `t`.

The orbit is assumed to be circular.
"""
orbital_velocity(sma::Unitful.Length, t::Unitful.Time) = 2π * sma / t |> u"km/s"
	
"""
    orbital_velocity(sma::Unitful.Length, t::Unitful.Time, e)
	
Approximate the average orbital velocity of a orbit with eccentricity `e` and period `t`.

``s_o = \\frac{2πa}{t} \\left( 1 - \\frac{e^2}{4} - \\frac{3e^4}{64} - \\frac{5e^6}{256} - \\frac{175e^8}{16384} \\right)``
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

``M_p = \\frac{4π^2a^3}{GT^2}``

where ``a`` is the semi-major axis and ``T`` is the orbital period.
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
	hill_sphere(m_parent::Unitful.Mass, m::Unitful.Mass, sma::Unitful.Length, e = 0)
	
Calculate the approximate Hill sphere radius at periapse of a body of mass `m` that orbits a body of mass `m_parent` at a distance of `sma` with orbital eccentricity `e`.

``h_r = a(1-e)\\sqrt[3]{\\frac{m}{3M}}``

where ``a`` is the semi-major axis, `m` is the mass of the smaller body and ``M`` is the mass of the central body.
"""
function hill_sphere(m_parent::Unitful.Mass, m::Unitful.Mass, sma::Unitful.Length, e = 0)
	sma * (1-e) * cbrt(m / 3m_parent) |> Unitful.unit(sma)
end

"""
	gravitational_binding_energy(m::Unitful.Mass, r::Unitful.Length)
	
Calculate the gravitational binding energy of a body with mass `m` and radius `r`.

``E_g = \\frac{3Gm^2}{5r}``
"""
gravitational_binding_energy(m::Unitful.Mass, r::Unitful.Length) = (3G*m^2) / 5r |>u"J"

"""
   roche_limit(r_primary::Unitful.Length, ρ_primary::Unitful.Density, ρ_satellite::Unitful.Density)
   
Rigid body approximation of the Roche limit for a body with radius `r_primary` and density `ρ_primary`, approached by a body with density `ρ_satellite`.

``r_r = r\\sqrt[2]{\frac{2ρ_p}{ρ_s}}``

where ``r`` and ``ρ_p`` are the radius and density of the central body, and ρ_s is the density of the orbiting body.
"""
roche_limit(r_primary::Unitful.Length, ρ_primary::Unitful.Density, ρ_satellite::Unitful.Density) = r_primary * cbrt(2ρ_primary / ρ_satellite) |> u"km"

"""
    gravity_tug_mass(force::Unitful.Force, separation::Unitful.Length, body_mass::Unitful.Mass)
	
Compute the mass required to apply `force` to a body of mass `body_mass` given their barycenters are at a distance of `separated`.
"""
gravity_tug_mass(force::Unitful.Force, separation::Unitful.Length, body_mass::Unitful.Mass) = (force * separation^2)/(G * body_mass) |> u"kg"

"""
    barycentric_distance(m_1::Unitful.Mass, m_2::Unitful.Mass, a::Unitful.Length)
	
Compute the distance from the barycenter of `m_1` to the barycenter of the `m_1`-`m_2` system, where the masses have an average separation of `a`.
"""
barycentric_distance(m1::Unitful.Mass, m2::Unitful.Mass, a::Unitful.Length) = (a * m2) / (m1 + m2) |> u"km"

"""
    tidal_acceleration(body_radius::Unitful.Length, parent_mass::Unitful.Mass, separation::Unitful.Length)
	
Approximate tidal acceleration felt by a body of radius `body_radius` at a distance of `separation` from a body of mass `parent_mass`.

Changes in tidal acceleration are associated with effects like heating, etc.

``a_t = \frac{2GMr}{s^3}``

where ``M`` is the mass of the central body, ``r`` is the radius of the body being affected and ``s`` is the distance between the two bodies.
"""
tidal_acceleration(body_radius::Unitful.Length, parent_mass::Unitful.Mass, separation::Unitful.Length) = 2body_radius * G * parent_mass / separation^3

function radial_orbit_time(m1::Unitful.Mass, m2::Unitful.Mass, x0::Unitful.Length, v0::Unitful.Velocity, x::Unitful.Length)
	μ = G * (m1 + m2)
	w = 1/x0 - v0^2/2μ
	
	return (asin(sqrt(w * x)) - sqrt(w * x * (1 - w * x))) / sqrt(2μ * w^3)
end

function radial_orbit_displacement(m1::Unitful.Mass, m2::Unitful.Mass, x0::Unitful.Length, v0::Unitful.Velocity, t::Unitful.Time)
	μ = G * (m1 + m2)
	w = 1/x0 - v0^2/2μ
	p = cbrt(9μ * t^2 / 2)
	
	return p - w*p^2 / 5 - 3w^2 * p^3 / 175 - 23w^3 * p^4/7875 - 1894w^4 * p^5 / 3931875 - 3293w^5 * p^6 / 21896875 - 2418092w^6 * p^7/62077640625
end

function radial_orbit_displacement2(m1::Unitful.Mass, m2::Unitful.Mass, x0::Unitful.Length, v0::Unitful.Velocity, t::Unitful.Time)
	μ = G * (m1 + m2)
	w = 1/x0 - v0^2/2μ
	p = cbrt(9μ * t^2 / 2)
	
	return p - w*p^2 / 5 - 3w^2 * p^3 / 175 - 23w^3 * p^4/7875 - 1894w^4 * p^5 / 3931875 - 3293w^5 * p^6 / 21896875 - 2418092w^6 * p^7/62077640625
end

function radial_orbit_displacement(m1::Real, m2::Real, x0::Real, v0::Real, t::Real)
	μ = ustrip(G) * (m1 + m2)
	w = 1/x0 - v0^2/2μ
	p = cbrt(9μ * t^2 / 2)
	
	return p - w*p^2 / 5 - 3w^2 * p^3 / 175 - 23w^3 * p^4/7875 - 1894w^4 * p^5 / 3931875 - 3293w^5 * p^6 / 21896875 - 2418092w^6 * p^7/62077640625
end

end
