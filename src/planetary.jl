module SfPlanetary

	using Unitful
	using UnitfulAstro
	using UnitfulAngles
	
	export Body, Orbit, Rotation

	const Degree{T} = Quantity{T, NoDims, typeof(u"°")}

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
		orbit::Union{Nothing, Orbit}
		Rotation::Union{Nothing, Rotation}
	end

	Orbit(parent::AbstractBody, semi_major_axis::Unitful.Length, eccentricity::Real, inclination::Degree) =
		Orbit(parent, semi_major_axis, eccentricity, nothing, inclination, nothing, nothing)
		
	Body(name::String, mass::Unitful.Mass, radius::Unitful.Length) =
		Body(name, mass, radius, radius, nothing, nothing)

	#####################################################################################
	
	import ..SfGravity: gravity, planetary_mass, planetary_radius, orbital_period, orbital_radius, escape_velocity, hill_sphere
	
	export gravity, planetary_mass, planetary_radius, orbital_period, orbital_radius, escape_velocity, hill_sphere

	gravity(body::Body) = gravity(body.mass, body.equatorial_radius)
	planetary_mass(body::Body) = body.mass
	planetary_radius(body::Body) = body.radius

	# kepler's third
	orbital_period(body::Body) = orbital_period(body.orbit.parent.mass, body.orbit.semi_major_axis)
	orbital_radius(body::Body) = body.orbit.semi_major_axis

	escape_velocity(body::Body) = escape_velocity(body.mass, body.equatorial_radius)

	hill_sphere(body::Body) = hill_sphere(body.orbit.parent.mass, body.mass, body.orbit.semi_major_axis, body.orbit.eccentricity)
	hill_sphere(orbit::Orbit, mass::Unitful.Mass) = hill_sphere(body.orbit.parent.mass, mass, body.orbit.semi_major_axis, body.orbit.eccentricity)

end