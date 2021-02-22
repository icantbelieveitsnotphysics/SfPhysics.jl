module SfLaser

using Unitful

export diffraction_limited_element_size, diffraction_limited_range, diffraction_limited_spot_diameter, diffraction_limited_wavelength
	attenuated_beam_power, approximate_breakdown_threshold, beam_intensity

function diffraction_limited_spot_diameter(range::Unitful.Length, element_diameter::Unitful.Length, wavelength::Unitful.Length)
	if range <= 0u"m"
		throw(DomainError(range, "Range must be greater than zero"))
	elseif element_diameter <= 0u"m"
		throw(DomainError(element_diameter, "Element diameter must be greater than zero"))
	elseif wavelength <= 0u"m"
		throw(DomainError(wavelength, "Wavelength must be greater than zero"))
	end

	if range < element_diameter
		@warn "Range less than element size, result will be inaccurate" range element_size
	end

	return 1.2range * wavelength / element_diameter |> u"m"
end

function diffraction_limited_element_size(range::Unitful.Length, spot_diameter::Unitful.Length, wavelength::Unitful.Length)
	if range <= 0u"m"
		throw(DomainError(range, "Range must be greater than zero"))
	elseif spot_diameter <= 0u"m"
		throw(DomainError(spot_diameter, "Spot diameter must be greater than zero"))
	elseif wavelength <= 0u"m"
		throw(DomainError(wavelength, "Wavelength must be greater than zero"))
	end

	element_diameter =  1.2range * wavelength / spot_diameter
	
	if range < element_diameter
		@warn "Range less than element size, result will be inaccurate" range element_diameter
	end
	
	return element_diameter |> u"m"
end

function diffraction_limited_range(element_diameter::Unitful.Length, spot_diameter::Unitful.Length, wavelength::Unitful.Length)
	if element_diameter <= 0u"m"
		throw(DomainError(element_diameter, "Element must be greater than zero"))
	elseif spot_diameter <= 0u"m"
		throw(DomainError(spot_diameter, "Spot diameter must be greater than zero"))
	elseif wavelength <= 0u"m"
		throw(DomainError(wavelength, "Wavelength must be greater than zero"))
	end

	range = (element_diameter * spot_diameter) / 1.2wavelength
	
	if range < element_diameter
		@warn "Range less than element size, result will be inaccurate" range element_diameter
	end
	
	return element_diameter |> u"m"
end

function diffraction_limited_wavelength(range::Unitful.Length, element_diameter::Unitful.Length, spot_diameter::Unitful.Length)
	if element_diameter <= 0u"m"
		throw(DomainError(element_diameter, "Element must be greater than zero"))
	elseif spot_diameter <= 0u"m"
		throw(DomainError(spot_diameter, "Spot diameter must be greater than zero"))
	elseif range <= 0u"m"
		throw(DomainError(range, "Range must be greater than zero"))
	end
	
	if range < element_diameter
		@warn "Range less than element size, result will be inaccurate" range element_diameter
	end
	
	return (element_diameter * spot_diameter) / range |> u"nm"
end

attenuated_beam_power(range::Unitful.Length, power::Unitful.Power, attenuation_length::Unitful.Length) = power * exp(-range / attenuation_length)

# http://panoptesv.com/SciFi/LaserDeathRay/Breakdown.html
function approximate_breakdown_threshold(wavelength::Unitful.Length, ambient_pressure::Unitful.Pressure)
	p = ustrip(ambient_pressure |> u"atm")
	l = ustrip(wavelength |> u"m")
	
	return ((0.31 / (l * p^0.6)^2) * u"W*cm^-2") |> u"W/m^2"
end

beam_intensity(power::Unitful.Power, spot_diameter::Unitful.Length) = power / (π * spot_diameter^2 / 4)

end