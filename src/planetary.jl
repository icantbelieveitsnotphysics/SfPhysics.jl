module SfPlanetary

	using Unitful
	using UnitfulAstro
	using UnitfulAngles
	
	export Body, Orbit, Rotation

	const Degree{T} = Quantity{T, NoDims, typeof(u"°")}
	@derived_dimension ThermalFlux Unitful.𝐌*Unitful.𝐓^-3

	abstract type AbstractBody end

	struct Orbit
		parent::AbstractBody
		semi_major_axis::Unitful.Length
		eccentricity::Real
		mean_anomaly::Union{Nothing, Degree}
		inclination::Degree # relative to equator of parent body
		ascending_node::Union{Nothing, Degree}
		periapsis::Union{Nothing, Degree}
	end

	struct Rotation
		moment_of_inertia::Union{Nothing, Real}
		rotation_period::Unitful.Time
		axial_tilt::Degree	
	end

	struct Body <: AbstractBody
		name::String
		mass::Unitful.Mass
		equatorial_radius::Unitful.Length
		polar_radius::Unitful.Length
		bond_albedo::Union{Nothing, Real}
		orbit::Union{Nothing, Orbit}
		rotation::Union{Nothing, Rotation}
	end

	Orbit(parent::AbstractBody, semi_major_axis::Unitful.Length, eccentricity::Real, inclination::Degree) =
		Orbit(parent, semi_major_axis, eccentricity, nothing, inclination, nothing, nothing)
		
	Body(name::String, mass::Unitful.Mass, radius::Unitful.Length, bond_albedo::Union{Nothing, Rotation} = nothing) =
		Body(name, mass, radius, radius, bond_albedo, nothing, nothing)

	#####################################################################################
	
	import ..SfGravity: gravity, planetary_mass, planetary_radius, orbital_period, orbital_radius, orbital_velocity, escape_velocity, hill_sphere
	import ..SfRelativity: relativistic_kinetic_energy
	import ..SfPhysics: kinetic_energy, spherical_cap_solid_angle
	
	import PhysicalConstants.CODATA2018: σ # Stefan-Boltzmann constant
	
	export gravity, planetary_mass, planetary_radius, orbital_period, orbital_radius, orbital_velocity, escape_velocity, hill_sphere,
		relativistic_kinetic_energy, kinetic_energy, stellar_irradiance, planetary_equilibrium_temperature

	"""
	gravity(body::Body)
	
Calculate the surface gravity of `body`.

# Example

```julia-repl
julia> gravity(SfSolarSystem.moon)
1.625143043265411976549419802968796018897730001404175282685973623864077610575064 m s^-2
```"""
	gravity(body::Body) = gravity(body.mass, body.equatorial_radius)
	planetary_mass(body::Body) = body.mass
	planetary_radius(body::Body) = body.equatorial_radius

	# kepler's third
	orbital_period(orbit::Orbit) = orbital_period(orbit.parent.mass, orbit.semi_major_axis)
	orbital_period(body::Body) = orbital_period(body.orbit)
	orbital_radius(orbit::Orbit) = orbit.semi_major_axis
	orbital_radius(body::Body) = orbital_radius(body.orbit)
	
	orbital_velocity(orbit::Orbit) = orbital_velocity(orbit.semi_major_axis, orbital_period(orbit), orbit.eccentricity)
	orbital_velocity(body::Body) = orbital_velocity(body.orbit)

	escape_velocity(body::Body) = escape_velocity(body.mass, body.equatorial_radius)

	hill_sphere(body::Body) = hill_sphere(body.orbit.parent.mass, body.mass, body.orbit.semi_major_axis, body.orbit.eccentricity)
	hill_sphere(orbit::Orbit, mass::Unitful.Mass) = hill_sphere(orbit.parent.mass, mass, orbit.semi_major_axis, orbit.eccentricity)

	kinetic_energy(mass::Unitful.Mass, orbit::Orbit) = kinetic_energy(mass, orbital_velocity(orbit))
	kinetic_energy(body::Body) = kinetic_energy(body.mass, body.orbit)
	
	"""
	stellar_luminosity(r_star::Unitful.Length, t_surface::Unitful.Temperature)
	
Approximate the luminosity of a star with radius `r_star` and surface temperature `t_surface`.
"""
	stellar_luminosity(r_star::Unitful.Length, t_surface::Unitful.Temperature) = 4π * r_star^2 * σ * t_star^4 |>u"W"
	
	"""
	stellar_irradiance = function(l_stellar::Unitful.Power, r_orbit::Unitful.Length, r_body::Unitful.Length)

Compute the proportion of a star's luminosity that falls upon a circular body of the given radius at the specified orbital distance.	
"""
	stellar_irradiance = function(l_stellar::Unitful.Power, r_orbit::Unitful.Length, r_body::Unitful.Length)
		Ω = spherical_cap_solid_angle(r_orbit, r_body) # steradians
		return l_stellar * Ω / 4π
	end
	
	planetary_equilibrium_temperature(irradiance::ThermalFlux, bond_albedo::Real) = ((irradiance * (1 - bond_albedo)) / 4σ)^.25 |> u"K"
	
	planetary_equilibrium_temperature(l_stellar::Unitful.Power, r_orbit::Unitful.Length, r_body::Unitful.Length, bond_albedo::Real) =
		planetary_equilibrium_temperature(solar_irradiance(s_l, r_orbit, r_body) / (π * r_body^2), bond_albedo) |> u"K"
		
	# temperature, exosphere altiutude, planetary mass, planetary radius, gas molecular mass
	# http://cococubed.asu.edu/code_pages/jeans_escape.shtml
	function jeans_escape_timescale(T::Unitful.Temperature, h::Unitful.Length, M::Unitful.Mass, R::Unitful.Length, m::Unitful.Mass)
       g = (G*M)/(R+h)^2
       H = (k_B*T)/(m*g) # scale height for gas
       v_peak=sqrt((2k_B*T)/m) # peak of maxwell-boltzmann distribution
       v_esc=sqrt((2G*M)/(R+h)) # escape velocity of planet at exosphere altitude
       λ =(v_esc/v_peak)^2
	   v_jeans = v_peak * ((1 + λ)*exp(-λ))/sqrt(4π)
       return H/v_jeans |> u"yr"
    end
end