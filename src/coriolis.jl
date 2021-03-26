module SfCoriolis

using Unitful, UnitfulAngles, Documenter
import LinearAlgebra: cross

import ..SfGeometry: radius, unit_x, unit_y, unit_z
import ..SfUnits: AngularVelocity, AngularMomentum, Angle, to_angle

export coriolis_acceleration, coriolis_force, tangential_velocity, centrifugal_acceleration, centrifugal_force, 
	angular_momentum, angular_velocity, radius
	
DocMeta.setdocmeta!(SfCoriolis, :DocTestSetup, :(using Unitful, UnitfulAstro, ..SfCoriolis, ..SfGeometry); recursive=true)

# this can't be a jldoctest example right now, as there's some regex/string comparison error that breaks it.
"""
    coriolis_acceleration(Ω::Vector, v::Vector)
	
The apparent acceleration due to coriolis effects for an object moving with velolcity
vector `v` in a rotating reference frame with angular velocity vector `Ω`.

# Example

```julia-repl
julia> coriolis_acceleration(0.001u"rad/s" * unit_y, 400u"m/s" * unit_x)
3-element Array{Quantity{Float64,𝐋  𝐓 ^-2,Unitful.FreeUnits{(m, s^-2),𝐋  𝐓 ^-2,nothing}},1}:
  0.0 m s^-2
  0.0 m s^-2
 -0.8 m s^-2
```
"""
coriolis_acceleration(Ω::Vector, v::Vector) = upreferred.(2cross(Ω, v))

"""
    coriolis_acceleration(Ω::Vector, v::Vector)
	
The strength of the coriolis force apparently acting on an object with mass `m` moving with velolcity
vector `v` in a rotating reference frame with angular velocity vector `Ω`.
"""
coriolis_force(m::Unitful.Mass, Ω::Vector, v::Vector) = mass * coriolis_force(Ω, v)

"""
    tangential_velocity(ω::AngularVelocity, r::Unitful.Length)
	
Calculate the tangential velocity for an object rotating about a point at a distance `r` with angular velocity `ω`.
"""
tangential_velocity(ω::AngularVelocity, r::Unitful.Length) = upreferred(ω * r)

"""
    tangential_velocity(ω::AngularVelocity, r::Unitful.Length, λ::Angle)
	
Calculate the tangential velocity of an object at latitude `λ` on a rotating sphere of radius `r` with angular velocity `ω`.
"""
tangential_velocity(ω::AngularVelocity, r::Unitful.Length, λ::Angle) = upreferred(ω * r * cos(λ))

"""
    centrifugal_acceleration(ω::AngularVelocity, r::Unitful.Length)
	
Calculate the centrifugal acceleration for an object rotating about a point at a distance `r` with angular velocity `ω`.
"""
centrifugal_acceleration(ω::AngularVelocity, r::Unitful.Length) = upreferred(ω^2 * r)

"""
    angular_velocity(ca::Unitful.Acceleration, r::Unitful.Length)
	
Calculate the angular velocity required to develop a centrifugal acceleration of `ca` at a radius of `r`.
"""
angular_velocity(ca::Unitful.Acceleration, r::Unitful.Length) = sqrt(ca / r) |> u"rad/s"

"""
    radius(ω::AngularVelocity, ca::Unitful.Acceleration)
	
Calculate the radius of rotation of a body with angular velocity `ω` to develop a centrifugal acceleration of `ca`.
"""
radius(ω::AngularVelocity, ca::Unitful.Acceleration) = upreferred(ca / ω^2)

"""
    centrifugal_acceleration(ω::AngularVelocity, r::Unitful.Length, λ::Angle)
	
Calculate the radial part of the centrifugal acceleration of object at latitude `λ` on a rotating sphere of radius `r` with angular velocity `ω`.
"""
centrifugal_acceleration(ω::AngularVelocity, r::Unitful.Length, λ::Angle) = upreferred(ω^2 * r * cos(λ))

"""
    centrifugal_force(m::Unitful.Mass, ω::AngularVelocity, r::Unitful.Length)
	
Calculate the centrifugal force acting on an object of mass `m` rotating about a point at a distance `r` with angular velocity `ω`.
"""
centrifugal_force(m::Unitful.Mass, ω::AngularVelocity, r::Unitful.Length) = upreferred(m * centrifugal_acceleration(ω, r))

"""
    angular_velocity(m::Unitful.Mass, cf::Unitful.Force, r::Unitful.Length)
	
Calculate the angular velocity required to develop a centrifugal force of `cf` on a body of mass `m` at a radius of `r`.
"""
angular_velocity(m::Unitful.Mass, cf::Unitful.Force, r::Unitful.Length) = angular_velocity(cf / m, r)

"""
    centrifugal_force(m::Unitful.Mass, ω::AngularVelocity, r::Unitful.Length, λ::Angle)
	
Calculate the radial part of the centrifugal force acting of object of mass `m` at latitude `λ` on a rotating sphere of radius `r` with angular velocity `ω`.
"""
centrifugal_force(m::Unitful.Mass, ω::AngularVelocity, r::Unitful.Length, λ::Angle) = upreferred(m * centrifugal_acceleration(ω, r, λ))

"""
    angular_momentum(m::Unitful.Mass, ω::AngularVelocity, r::Unitful.Length)
	
Calculate the angular momentum of an object of mass `m` rotating about a point at a distance `r` with angular velocity `ω`.
"""
angular_momentum(m::Unitful.Mass, ω::AngularVelocity, r::Unitful.Length) = upreferred(m * r^2 * ω)

"""
	angular_velocity(m::Unitful.Mass, L::AngularMomentum, r::Unitful.Length)
	
Calculate the angular velocity of an object with mass `m` and angular momentum `L` rotating about a point at a distance of `r`.
"""
angular_velocity(m::Unitful.Mass, L::AngularMomentum, r::Unitful.Length) = L / (m * r^2) |> u"rad/s"

end
