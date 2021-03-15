module SfElectro

using Unitful

import PhysicalConstants.CODATA2018: μ_0

export energy_density, field_strength, dipole_distance

@derived_dimension MagneticDipole Unitful.𝐈*Unitful.𝐋^2

"""
    energy_density(b::Unitful.BField)
	
Energy density of a uniform magnetic field
"""
energy_density(b::Unitful.BField) = b^2 / 2μ_0 |> u"J/m^3"

"""
    field_strength(ed::Unitful.Quantity)
	
Given a uniform magnetic field of energy density `ed`, compute the required field strength.
"""
field_strength(ed::Unitful.Quantity) = sqrt(ed * 2μ_0) |> u"T"

"""
    lorentz_force(i::Unitful.Current, l::Unitful.Length, b::Unitful.BField)
	
Force exerted on a charge carrying wire of length `l` with current `i` in a perpendicular magnetic field of strength `b`.
"""
lorentz_force(i::Unitful.Current, l::Unitful.Length, b::Unitful.BField) = i * l * b |>u"N"

"""
    field_strength(m::DipoleStrength, r::Unitful.Length)
	
Field strength on the equatorial plane of a magnetic dipole `m` at a distance `r` from the magnetic axis.

This assumes that the scale of the dipole is small compared to `r`.
"""
field_strength(m::MagneticDipole, r::Unitful.Length) = (2m * μ_0) / (4π * r^3) |> u"T"

"""
    dipole_distance(m::MagneticDipole, b::Unitful.BField)
	
Compute the distance away from the axis of magnetic dipole `m` on its equatorial plane that a field of strength `b` would be found.
"""
dipole_distance(m::MagneticDipole, b::Unitful.BField) = cbrt((2m * μ_0) / (4π * b)) |> u"m"

end
