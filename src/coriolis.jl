module SfCoriolis

using Unitful, UnitfulAngles
import LinearAlgebra: cross

export coriolis

const AngularVelocity{T} = Quantity{T, Unitful.𝐓 ^-1, typeof(u"rad/s")}

coriolis_force(rotation::Vector{AngularVelocity}, velocity::Vector{Unitful.Velocity}) = 2cross(rotation, velocity)

coriolis_force(mass::Unitful.Mass, rotation::Vector{AngularVelocity}, velocity::Vector{Unitful.Velocity}) = mass * coriolis_force(rotation, velocity)

#centrifugal_force(mass::

end