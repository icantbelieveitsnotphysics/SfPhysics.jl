module SfPlanetary

using Unitful, UnitfulAstro, UnitfulAngles
using ..SfGeometry, ..SfThermo

export Body, Orbit, Rotation, Star
export satellites

import ..SfUnits: Angle, to_angle

@derived_dimension ThermalFlux Unitful.𝐌*Unitful.𝐓^-3

abstract type AbstractBody end

const satellites_of = Dict{AbstractBody, Vector{AbstractBody}}()

"""
    satellites(body::AbstractBody)
	
Return a list of sattelites associated with `body`, which may be empty.
"""
satellites(body::AbstractBody) = haskey(satellites_of, body) ? satellites_of[body] : Vector{AbstractBody}()

"""
    add_satellite(primary::AbstractBody, satellite::AbstractBody)
	
Associate `satellite` with `primary`. Should not usually be called by user code.
"""
function add_satellite(primary::AbstractBody, satellite::AbstractBody)
	if !haskey(satellites_of, primary)
		satellites_of[primary] = [ satellite ]
	else
		push!(satellites_of[primary], satellite)
	end
end

struct Orbit
	parent::AbstractBody
	semi_major_axis::Unitful.Length
	eccentricity::Real
	inclination::Angle # relative to equator of parent body
	ascending_node::Union{Missing, Angle}
	periapsis::Union{Missing, Angle}
end

Orbit(parent::AbstractBody, 
		semi_major_axis::Unitful.Length, eccentricity::Real, 
		inclination::Real, ascending_node::Union{Missing, Real} = missing, periapsis::Union{Missing, Real} = missing) =
	Orbit(parent, semi_major_axis, eccentricity, to_angle(inclination), to_angle(ascending_node), to_angle(periapsis))
	
# needed to avoid some awkward ambiguous method resolution issues :(
Orbit(parent::AbstractBody, semi_major_axis::Unitful.Length, eccentricity::Real, inclination::Angle) =
	Orbit(parent, semi_major_axis, eccentricity, to_angle(inclination), missing, missing)

"""
    Rotation(moment_of_inertia::Union{Nothing, Real}, rotation_period::Unitful.Time, axial_tilt::Angle)
	
Describe the rotation of a body about its axis, with an optional `moment_of_inertia` and an `axial_tilt` relative to either its orbit or the local ecliptic depending on context.
"""
struct Rotation
	moment_of_inertia::Union{Missing, Real}
	rotation_period::Unitful.Time
	axial_tilt::Angle	
end

struct Star <: AbstractBody
	name::String
	mass::Unitful.Mass
	shape::Shape
	surface_temperature::Unitful.Temperature
	absolute_magnitude::Real
	spectral_class::String
	
	orbit::Union{Nothing, Orbit}
	rotation::Union{Nothing, Rotation}
	
	function Star(name::String, mass::Unitful.Mass, shape::Shape, temp::Unitful.Temperature, amag::Real, class::String, orbit::Union{Nothing, Orbit} = nothing, rotation::Union{Nothing, Rotation} = nothing)
		b = new(name, mass, shape, temp, amag, class, orbit, rotation)
		
		if orbit != nothing
			add_satellite(orbit.parent, b)
		end
		
		return b
	end
end

"""
    Star(name::String, mass::Unitful.Mass, radius::Unitful.Length, temp::Unitful.Temperature, amag::Real, class::String)
	
Construct a spherical, non-rotating, solitary star with the given characteristics.
"""
Star(name::String, mass::Unitful.Mass, radius::Unitful.Length, temp::Unitful.Temperature, amag::Real, class::String) = Star(name, mass, Sphere(radius), temp, amag, class)

struct Body <: AbstractBody
	name::String
	mass::Unitful.Mass
	shape::Shape
	bond_albedo::Union{Missing, Real}
	geometric_albedo::Union{Missing, Real}
	orbit::Union{Nothing, Orbit}
	rotation::Union{Nothing, Rotation}
	
"""
    Body(name::String, mass::Unitful.Mass, shape::Shape, bond_albedo::Union{Missing, Real}, geometric_albedo::Union{Missing, Real}, orbit::Union{Nothing, Orbit}, rotation::Union{Nothing, Rotation})
	
Construct a new body with the given characteristics, with an optional `orbit` defining its relationship with a parent body and optional `rotation` defining the nature of its day.
"""
	function Body(name::String, mass::Unitful.Mass, shape::Shape, bond_albedo::Union{Missing, Real}, geometric_albedo::Union{Missing, Real}, orbit::Union{Nothing, Orbit}, rotation::Union{Nothing, Rotation})
		b = new(name, mass, shape, bond_albedo, geometric_albedo, orbit, rotation)
		
		if orbit != nothing
			add_satellite(orbit.parent, b)
		end
		
		return b
	end
end

"""
    Body(name::String, mass::Unitful.Mass, radius::Unitful.Length, bond_albedo::Union{Nothing, Real} = nothing)
	
Construct a spherical, non-rotating, solitary astronomical body with the given characteristics.
"""	
Body(name::String, mass::Unitful.Mass, radius::Unitful.Length, bond_albedo::Union{Missing, Real} = missing, geometric_albedo::Union{Missing, Real} = missing) =
	Body(name, mass, Sphere(radius), bond_albedo, geometric_albedo, nothing, nothing)
	
Rotation(moment_of_inertia::Union{Missing, Real}, rotation_period::Unitful.Time, axial_tilt::Real) =
	Rotation(moment_of_inertia, rotation_period, to_angle(axial_tilt))

#####################################################################################

import ..SfGravity: gravity, planetary_mass, planetary_radius, orbital_period, orbital_radius, orbital_velocity, escape_velocity, 
	hill_sphere, gravitational_binding_energy, roche_limit, barycentric_distance
import ..SfRelativity: relativistic_kinetic_energy
import ..SfPhysics: kinetic_energy
import ..SfMatter: density, mass
import ..SfGeometry: spherical_cap_solid_angle, volume, radius, equatorial_radius, polar_radius, area, cross_sectional_area
import ..SfAstronomy: absolute_magnitude, apparent_magnitude, diffuse_sphere_q

import PhysicalConstants.CODATA2018: σ, G, k_B # σ = Stefan-Boltzmann constant, k_B Boltzmann constant

export gravity, orbital_period, orbital_radius, orbital_velocity, escape_velocity, hill_sphere,
	relativistic_kinetic_energy, kinetic_energy, stellar_luminosity, stellar_irradiance, planetary_equilibrium_temperature,
	jeans_escape_timescale, jeans_parameter, gravitational_binding_energy, roche_limit, volume, density, radius, equatorial_radius,
	area, cross_sectional_area, absolute_magnitude, apparent_magnitude, star, stellar_distance

"""
	gravity(body::AbstractBody)
	
Calculate the surface gravity of `body`.

# Example

```julia-repl
julia> gravity(SfSolarSystem.moon)
1.625143043265411976549419802968796018897730001404175282685973623864077610575064 m s^-2
```"""
gravity(body::AbstractBody) = gravity(body.mass, equatorial_radius(body))
planetary_mass(body::AbstractBody) = body.mass
planetary_radius(body::AbstractBody) = equatorial_radius(body)

# kepler's third
orbital_period(orbit::Orbit) = orbital_period(orbit.parent.mass, orbit.semi_major_axis)
orbital_period(body::AbstractBody) = orbital_period(body.orbit)
orbital_radius(orbit::Orbit) = orbit.semi_major_axis
orbital_radius(body::AbstractBody) = orbital_radius(body.orbit)

orbital_velocity(orbit::Orbit) = orbital_velocity(orbit.semi_major_axis, orbital_period(orbit), orbit.eccentricity)
orbital_velocity(body::AbstractBody) = orbital_velocity(body.orbit)

"""
    escape_velocity(body::AbstractBody)
	
Calculate the escape velocity from the surface of `body`.

If `body` is spheroidal, the escaping object is taken to be on its equator, otherwise the average radius of `body` is used.
"""
escape_velocity(body::AbstractBody) = escape_velocity(body.mass, radius(body))

"""
    escape_velocity(orbit::Orbit)
	
Calculate the escape velocity at the given `orbit`, based on the mass of the primary and the semimajor axis.
"""
escape_velocity(orbit::Orbit) = escape_velocity(orbit.parent.mass, orbit.semi_major_axis)

"""
    hill_sphere(body::AbstractBody)
	
Calculate the radius of the Hill sphere of `body`, using its current orbit and parent.
"""
hill_sphere(body::AbstractBody) = hill_sphere(body.orbit.parent.mass, body.mass, body.orbit.semi_major_axis, body.orbit.eccentricity)

"""
    hill_sphere(orbit::Orbit, mass::Unitful.Mass)
	
Calculate the radius of the Hill sphere of a body of `mass` in `orbit`, using the parent body described by that orbit.
"""
hill_sphere(orbit::Orbit, mass::Unitful.Mass) = hill_sphere(orbit.parent.mass, mass, orbit.semi_major_axis, orbit.eccentricity)

"""
    hill_sphere(parent::AbstractBody, satellite::AbstractBody, sma::Unitful.Length, ecc::Real = 1)
	
Calculate the radius of the Hill sphere of `satellite` as if it were orbiting `parent` with a semi-major axis of `sma` and eccentricity `ecc`.
"""
hill_sphere(parent::AbstractBody, satellite::AbstractBody, sma::Unitful.Length, ecc::Real = 1) = hill_sphere(mass(parent), mass(satellite), sma, ecc)

"""
    hill_sphere(parent::AbstractBody, satellite::AbstractBody, orbit::Orbit)
	
Calculate the radius of the Hill sphere of `satellite` as if it were orbiting `parent` in orbit `orbit`.
"""
hill_sphere(parent::AbstractBody, satellite::AbstractBody, orbit::Orbit) = hill_sphere(mass(parent), mass(satellite), orbit.semi_major_axis, orbit.eccentricity)

"""
    gravitational_binding_energy(body::AbstractBody)
	
Calculate the gravitational binding energy of `body`.	
"""
gravitational_binding_energy(body::AbstractBody) = gravitational_binding_energy(body.mass, radius(body))

"""
    kinetic_energy(mass::Unitful.Mass, orbit::Orbit)
	
Compute the average kinetic energy of a body of mass `mass` in orbit `orbit`.
"""
kinetic_energy(mass::Unitful.Mass, orbit::Orbit) = kinetic_energy(mass, orbital_velocity(orbit))

"""
    kinetic_energy(body::AbstractBody)
	
Compute the average kinetic energy of `body` in its current orbit.
"""
kinetic_energy(body::AbstractBody) = kinetic_energy(body.mass, body.orbit)

"""
    roche_limit(primary::AbstractBody, satellite::AbstractBody)
	
Calculate the Roche limit of `satellite` as it approaches `primary`.
"""
roche_limit(primary::AbstractBody, satellite::AbstractBody) = roche_limit(equatorial_radius(primary), density(primary), density(satellite))

"""
    roche_limit(body::AbstractBody)
	
Calculate the Roche limit of `body` as it approache the object it orbits.

If `body` does not orbit anything, an error will be raised.
"""
roche_limit(body::AbstractBody) = roche_limit(equatorial_radius(body.orbit.parent), density(body.orbit.parent), density(body))

"""
    barycentric_distance(body::AbstractBody)
	
Find the offset from the barycentre of `body` and its parent to the centre of the parent.

If `body` does not orbit anything, an error will be raised.
"""
barycentric_distance(body::AbstractBody) = barycentric_distance(mass(body), mass(body.orbit.parent), orbital_radius(body))

"""
    volume(body::AbstractBody)
	
Calculate the volume of `body`.
"""
volume(body::AbstractBody) = volume(body.shape) |> u"km^3"

"""
    area(body::AbstractBody)
	
Calculate the surface area of `body`.

Some shapes, such as triaxial ellipsoids, may not be supported yet.
"""
area(body::AbstractBody) = area(body.shape) |> u"km^2"

"""
    cross_sectional_area(body::AbstractBody)
	
Calculate the approximate cross sectional area of `body`, where possible.

Results are well defined for spherical and spheroidal bodies, but may throw an error for other shapes such as triaxial ellipsoids.
"""
cross_sectional_area(body::AbstractBody) = cross_sectional_area(body.shape) |> u"km^2"

"""
    radius(body::AbstractBody)
	
Returns the equatorial radius of `body` when it is spheroidal, or its average radius if it is ellipsoidal.
"""
radius(body::AbstractBody) = is_triaxial(body.shape) ? radius(body.shape) : equatorial_radius(body.shape)

"""
    equatorial_radius(body::AbstractBody)
	
Returns the equatorial radius of `body` when it is spheroidal, or an error for bodies without a well defined equatorial radius.
"""
equatorial_radius(body::AbstractBody) = equatorial_radius(body.shape)

"""
    polar_radius(body::AbstractBody)
	
Returns the polar radius of `body` when it is spheroidal, or an error for bodies without a well defined polar radius.
"""
polar_radius(body::AbstractBody) = polar_radius(body.shape)

"""
    radius(orbit::Orbit)
	
Returns the semi-major axis of `orbit`.
"""
radius(orbit::Orbit) = orbit.semi_major_axis

"""
    mass(body::AbstractBody)
	
Return the mass of `body`.
"""
mass(body::AbstractBody) = body.mass

"""
    star(body::AbstractBody)
	
Recuse upwards through the hierarchy of orbits to find the star that `body` ultimately orbits. If there is no parent star, `nothing` is returned.
"""
function star(body::AbstractBody)
	if typeof(body) == Star
		return star
	elseif body.orbit == nothing
		return nothing
	else
		return star(body.orbit.parent)
	end
end

"""
	stellar_distance(body::AbstractBody)
	
Returns the approximate average distance of `body` to its parent star, if any. If there is no parent star, an error is raised. If `body` is a star, 0 is returned.

This returns the average radius of the top-level orbit, and so is only a loose approximation for eccentric.
"""
function stellar_distance(body::AbstractBody)
	if typeof(body) == Star
		return 0u"AU"
	elseif body.orbit == nothing
		error("Body does not ultimately orbit a star")
	elseif typeof(body.orbit.parent) == Star
		return radius(body.orbit)
	else
		return stellar_distance(body.orbit.parent)
	end
end
	
"""
	stellar_luminosity(r_star::Unitful.Length, t_surface::Unitful.Temperature, r_star_p = r_star)
	
Approximate the luminosity of a star with radius `r_star` and surface temperature `t_surface`.

Optionally, a polar radius may be given, otherwise the star is assumed to have a circular cross-section.
"""
stellar_luminosity(r_star::Unitful.Length, t_surface::Unitful.Temperature, r_star_p = r_star) = 4π * r_star * r_star_p * σ * t_surface^4 |>u"W"

"""
    stellar_luminosity(star::Star)
	
Approximate the luminosity of the given `star`.
"""
stellar_luminosity(star::Star) = stellar_luminosity(equatorial_radius(star), star.surface_temperature, polar_radius(star))
	
"""
	stellar_irradiance = function(l_stellar::Unitful.Power, r_orbit::Unitful.Length, r_body::Unitful.Length)

Compute the proportion of a star's luminosity that falls upon a circular body of the given radius at the specified orbital distance.	
"""
function stellar_irradiance(l_stellar::Unitful.Power, r_orbit::Unitful.Length, r_body::Unitful.Length)
	Ω = spherical_cap_solid_angle(r_orbit, r_body) # steradians
	return l_stellar * Ω / 4π
end

"""
	stellar_irradiance(body::AbstractBody)
	
Approximate the irradiance delivered to `body` by the star it ultimately orbits.

This function will recursively ascend its tree of parent bodies until a star is found. If the ultimate parent is not
a star, an error will be raised. The distance from `star` to the topmost parent is used to compute the irradiance, 
as a proxy for average orbital distance.

This method will deadlock if used on a binary or multiple system at present.
"""
function stellar_irradiance(body::AbstractBody)
	s = star(body)
		
	if s == nothing
		error("Body does not ultimately orbit a star")
	end
	
	return stellar_irradiance(stellar_luminosity(s), orbit.semi_major_axis, radius(body))
end

"""
    planetary_equilibrium_temperature(irradiance::ThermalFlux, bond_albedo::Real)
	
Compute the equilibrium temperature for a body with `bond_albedo` subject to stellar `irradiance`.
"""
planetary_equilibrium_temperature(irradiance::ThermalFlux, bond_albedo::Real) = ((irradiance * (1 - bond_albedo)) / 4σ)^.25 |> u"K"

"""
    planetary_equilibrium_temperature(l_stellar::Unitful.Power, r_orbit::Unitful.Length, r_body::Unitful.Length, bond_albedo::Real)
	
Approximate the planetary equilibrium temperature of a body of radius `r_body` and `bond_albedo`, orbiting a star of luminosity `l_stellar` at a distance of `r_orbit`.
"""
planetary_equilibrium_temperature(l_stellar::Unitful.Power, r_orbit::Unitful.Length, r_body::Unitful.Length, bond_albedo::Real) =
	planetary_equilibrium_temperature(stellar_irradiance(s_l, r_orbit, r_body) / (π * r_body^2), bond_albedo) |> u"K"

"""
    planetary_equilibrium_temperature(body::Body)
	
Approximate the planetary equilibrium temperature of `body`, using the irradiance of the star it ultimately orbits.

If the body does not ultimately orbit a star, an error will be raised. Binary and multiple star systems will result in deadlocks.
"""
planetary_equilibrium_temperature(body::Body) = 
	planetary_equilibrium_temperature(stellar_irradiance(body) / cross_sectional_area(body), body.bond_albedo) |> u"K"
	
# temperature, exosphere altiutude, planetary mass, planetary radius, gas molecular mass
# http://cococubed.asu.edu/code_pages/jeans_escape.shtml
function jeans_escape_timescale(T::Unitful.Temperature, h::Unitful.Length, M::Unitful.Mass, R::Unitful.Length, m::Unitful.Mass)
   g = (G*M)/(R+h)^2 # gravity at exobase
   H = (k_B*T)/(m*g) # scale height for gas
   v_peak=maxwell_boltzmann_peak_speed(m, T) # peak of maxwell-boltzmann distribution
   v_esc=sqrt((2G*M)/(R+h)) # escape velocity of planet at exosphere altitude
   λ =(v_esc/v_peak)^2
   v_jeans = v_peak * ((1 + λ)*exp(-λ))/sqrt(4π)
   return H/v_jeans |> u"yr"
end

jeans_escape_timescale(T::Unitful.Temperature, h::Unitful.Length, body::AbstractBody, m::Unitful.Mass) = 
	jeans_escape_timescale(T, h, mass(body), radius(body), m)

# https://arxiv.org/ftp/arxiv/papers/1009/1009.5110.pdf
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008GL036513
jeans_parameter(m_planet::Unitful.Mass, m_molecule::Unitful.Mass, t_exosphere::Unitful.Temperature, r_exosphere::Unitful.Length) = 
	(G * m_planet * m_molecule) / (k_B * t_exosphere * r_exosphere) |> u"m/m"

jeans_parameter(body::AbstractBody, m_molecule::Unitful.Mass, t_exosphere::Unitful.Temperature, r_exosphere::Unitful.Length) = 
	jeans_parameter(mass(body), m_molecule, t_exosphere, r_exosphere)
	
"""
    absolute_magnitude(body::Body)
	
Absolute magnitude of `body`, assuming it is a diffusely reflecting sphere with a known geometric albedo.
"""
absolute_magnitude(body::Body) = absolute_magnitude(radius(body), body.geometric_albedo)

"""
    apparent_magnitude(observer::AbstractBody, observed::Body, separation::Unitful.Length, phase::Angle)
	
Apparent magnitude of `observed` as seen from `observer` at a distance of `separation` with `phase` angle.

`observed` is assumed to be a diffusely reflecting sphere. Values will be incorrect for non-spherical bodies.
"""
apparent_magnitude(observer::AbstractBody, observed::Body, separation::Unitful.Length, phase::Angle) =
	apparent_magnitude(absolute_magnitude(observed), stellar_distance(observed), separation, stellar_distance(observer), diffuse_sphere_q(phase))

end
