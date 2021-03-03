using Unitful
using .SfElectro

import .SfElectro.Tesla

export cratering_depth, cratering_strength, cratering_volume, jet_penetrator, newtonian_penetrator, coilgun_velocity

function cratering_strength(yield_strength::Unitful.Pressure)
    3u"J/m^3" * Unitful.ustrip(u"Pa", yield_strength)
end

function cratering_volume(e::Unitful.Energy, yield_strength::Unitful.Pressure)
    e / cratering_strength(yield_strength)
end

function cratering_depth(e::Unitful.Energy, yield_strength::Unitful.Pressure)
    r = cratering_volume(e, yield_strength)
    # hemisphere volume = (2/3)πr^3

    cbrt(3r/2π)
end

"""
    jet_penetrator(l::Unitful.Length, ρj::Unitful.Density, ρt::Unitful.Density)
	
Approximation for penetration depth of hypervelocity jet of density `ρj` and length `l` into solid material of density `ρt`.
"""
jet_penetrator(l::Unitful.Length, ρj::Unitful.Density, ρt::Unitful.Density) = l * sqrt(ρj / ρt)

"""
    newtonian_penetrator(l::Unitful.Length, ρj::Unitful.Density, ρt::Unitful.Density)
	
Approximation for penetration depth of high velocity impactor of density `ρj` and length `l` into solid material of density `ρt`.
"""
newtonian_penetrator(l::Unitful.Length, ρj::Unitful.Density, ρt::Unitful.Density) = l * (ρj / ρt)

"""
    coilgun_velocity(m::Unitful.Mass, r::Unitful.Length, l::Unitful.Length, b::Tesla)
	
Use Luke Campbell's coilgun approximation to compute muzzle velocity of a cylindrical projectile of mass `m` and radius `r` firing down a barrel of length `l` in a unitform magnetic field of strength `b`.
"""
function coilgun_velocity(m::Unitful.Mass, r::Unitful.Length, l::Unitful.Length, b::Tesla)
	ed = energy_density(b)
	v = π * r^2 * l
	e = ed * v
	
	return sqrt(2e / m) |> u"km/s"
end