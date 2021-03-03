module SfElectro

using Unitful

import PhysicalConstants.CODATA2018: μ_0

export energy_density

Tesla = typeof(0u"T")

"""
    energy_density(b::Tesla)
	
Energy density of a magnetic field
"""
energy_density(b::Tesla) = b^2 / 2μ_0 |> u"J/m^3"

end