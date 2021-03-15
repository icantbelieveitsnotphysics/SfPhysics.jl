module SfPlanetary

using Unitful, UnitfulAstro, UnitfulAngles
using ..SfGeometry

export Body, Orbit, Rotation, Star
export satellites

import ..SfPhysics: Angle, to_angle

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
	ascending_node::Union{Nothing, Angle}
	periapsis::Union{Nothing, Angle}
end

"""
    Rotation(moment_of_inertia::Union{Nothing, Real}, rotation_period::Unitful.Time, axial_tilt::Angle)
	
Describe the rotation of a body about its axis, with an optional `moment_of_inertia` and an `axial_tilt` relative to either its orbit or the local ecliptic depending on context.
"""
struct Rotation
	moment_of_inertia::Union{Nothing, Real}
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
	bond_albedo::Union{Nothing, Real}
	orbit::Union{Nothing, Orbit}
	rotation::Union{Nothing, Rotation}
	
"""
    Body(name::String, mass::Unitful.Mass, shape::Shape, bond_albedo::Union{Nothing, Real}, orbit::Union{Nothing, Orbit}, rotation::Union{Nothing, Rotation})
	
Construct a new body with the given characteristics, with an optional `orbit` defining its relationship with a parent body and optional `rotation` defining the nature of its day.
"""
	function Body(name::String, mass::Unitful.Mass, shape::Shape, bond_albedo::Union{Nothing, Real}, orbit::Union{Nothing, Orbit}, rotation::Union{Nothing, Rotation})
		b = new(name, mass, shape, bond_albedo, orbit, rotation)
		
		if orbit != nothing
			add_satellite(orbit.parent, b)
		end
		
		return b
	end
end

Orbit(parent::AbstractBody, 
		semi_major_axis::Unitful.Length, eccentricity::Real, 
		inclination::Real, ascending_node::Union{Nothing, Real} = nothing, periapsis::Union{Nothing, Real} = nothing) =
	Orbit(parent, semi_major_axis, eccentricity, to_angle(inclination), to_angle(ascending_node), to_angle(periapsis))
	
# needed to avoid some awkward ambiguous method resolution issues :(
Orbit(parent::AbstractBody, semi_major_axis::Unitful.Length, eccentricity::Real, inclination::Angle) =
	Orbit(parent, semi_major_axis, eccentricity, to_angle(inclination), nothing, nothing)

"""
    Body(name::String, mass::Unitful.Mass, radius::Unitful.Length, bond_albedo::Union{Nothing, Real} = nothing)
	
Construct a spherical, non-rotating, solitary astronomical body with the given characteristics.
"""	
Body(name::String, mass::Unitful.Mass, radius::Unitful.Length, bond_albedo::Union{Nothing, Real} = nothing) =
	Body(name, mass, Sphere(radius), bond_albedo, nothing, nothing)
	
Rotation(moment_of_inertia::Union{Nothing, Real}, rotation_period::Unitful.Time, axial_tilt::Real) =
	Rotation(moment_of_inertia, rotation_period, to_angle(axial_tilt))

#####################################################################################

import ..SfGravity: gravity, planetary_mass, planetary_radius, orbital_period, orbital_radius, orbital_velocity, escape_velocity, 
	hill_sphere, gravitational_binding_energy, roche_limit
import ..SfRelativity: relativistic_kinetic_energy
import ..SfPhysics: kinetic_energy
import ..SfMatter: density, mass
import ..SfGeometry: spherical_cap_solid_angle, volume, radius, equatorial_radius, area

import PhysicalConstants.CODATA2018: σ, G, k_B # σ = Stefan-Boltzmann constant, k_B Boltzmann constant

export gravity, planetary_mass, planetary_radius, orbital_period, orbital_radius, orbital_velocity, escape_velocity, hill_sphere,
	relativistic_kinetic_energy, kinetic_energy, stellar_luminosity, stellar_irradiance, planetary_equilibrium_temperature,
	jeans_escape_timescale, jeans_parameter, gravitational_binding_energy, roche_limit, volume, density, radius, equatorial_radius

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
planetary_radius(body::AbstractBody) = body.equatorial_radius

# kepler's third
orbital_period(orbit::Orbit) = orbital_period(orbit.parent.mass, orbit.semi_major_axis)
orbital_period(body::AbstractBody) = orbital_period(body.orbit)
orbital_radius(orbit::Orbit) = orbit.semi_major_axis
orbital_radius(body::AbstractBody) = orbital_radius(body.orbit)

orbital_velocity(orbit::Orbit) = orbital_velocity(orbit.semi_major_axis, orbital_period(orbit), orbit.eccentricity)
orbital_velocity(body::AbstractBody) = orbital_velocity(body.orbit)

escape_velocity(body::AbstractBody) = escape_velocity(body.mass, body.equatorial_radius)

hill_sphere(body::AbstractBody) = hill_sphere(body.orbit.parent.mass, body.mass, body.orbit.semi_major_axis, body.orbit.eccentricity)
hill_sphere(orbit::Orbit, mass::Unitful.Mass) = hill_sphere(orbit.parent.mass, mass, orbit.semi_major_axis, orbit.eccentricity)

gravitational_binding_energy(body::AbstractBody) = gravitational_binding_energy(body.mass, body.equatorial_radius)

kinetic_energy(mass::Unitful.Mass, orbit::Orbit) = kinetic_energy(mass, orbital_velocity(orbit))
kinetic_energy(body::AbstractBody) = kinetic_energy(body.mass, body.orbit)

roche_limit(primary::AbstractBody, satellite::AbstractBody) = roche_limit(primary.equatorial_radius, density(primary), density(satellite))

volume(body::AbstractBody) = volume(body.shape) |> u"km^3"

area(body::AbstractBody) = area(body.shape) |> u"km^2"

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
    radius(orbit::Orbit)
	
Returns the semi-major axis of `orbit`.
"""
radius(orbit::Orbit) = orbit.semi_major_axis

mass(body::AbstractBody) = body.mass
	
"""
	stellar_luminosity(r_star::Unitful.Length, t_surface::Unitful.Temperature)
	
Approximate the luminosity of a star with radius `r_star` and surface temperature `t_surface`.
"""
stellar_luminosity(r_star::Unitful.Length, t_surface::Unitful.Temperature) = 4π * r_star^2 * σ * t_surface^4 |>u"W"

"""
    stellar_luminosity(star::Star)
	
Approximate the luminosity of the given `star`.
"""
stellar_luminosity(star::Star) = stellar_luminosity(radius(star), star.surface_temperature)
	
"""
	stellar_irradiance = function(l_stellar::Unitful.Power, r_orbit::Unitful.Length, r_body::Unitful.Length)

Compute the proportion of a star's luminosity that falls upon a circular body of the given radius at the specified orbital distance.	
"""
function stellar_irradiance(l_stellar::Unitful.Power, r_orbit::Unitful.Length, r_body::Unitful.Length)
	Ω = spherical_cap_solid_angle(r_orbit, r_body) # steradians
	return l_stellar * Ω / 4π
end

"""
	stellar_irradiance(star::Star, body::AbstractBody)
	
Approximate the irradiance delivered by `star` to `body`.

If `body` does not directly orbit `star`, this function will recursively ascend its tree of parents until the star is found, 
or it is discovered not to belong to `star`'s system which will raise an error. The distance from `star` to the topmost parent
is used to compute the irradiance, as a proxy for average orbital distance.

This method will deadlock if used on a binary system at present.
"""
function stellar_irradiance(star::Star, body::AbstractBody)
	orbit = body.orbit
	
	while orbit.parent != star
		orbit = orbit.parent.orbit
		
		if orbit == nothing
			error("Body does not ultimately orbit given star")
		end
	end
	
	return stellar_irradiance(stellar_luminosity(star), orbit.semi_major_axis, radius(body))
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

# https://arxiv.org/ftp/arxiv/papers/1009/1009.5110.pdf
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008GL036513
jeans_parameter(m_planet::Unitful.Mass, m_molecule::Unitful.Mass, t_exosphere::Unitful.Temperature, r_exosphere::Unitful.Length) = 
	(G * m_planet * m_molecule) / (k_B * t_exosphere * r_exosphere) |> u"m/m"

end
