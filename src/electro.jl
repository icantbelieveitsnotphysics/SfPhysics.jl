module SfElectro

using Unitful

import PhysicalConstants.CODATA2018: μ_0

export energy_density, field_strength

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

end
