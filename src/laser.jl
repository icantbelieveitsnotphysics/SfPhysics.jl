module SfLaser

using Unitful

export diffraction_limited_element_size, diffraction_limited_spot_diameter, diffraction_limited_range, diffraction_limited_wavelength,
	attenuated_beam_power, approximate_breakdown_threshold, beam_intensity
	
"""
    diffraction_limited_spot_diameter(r::Unitful.Length, d_element::Unitful.Length, λ::Unitful.Length)
	
Approximate minimum spot size for a diffraction-limited gaussian laser beam of wavelength `λ` with
initial diameter `d_element` at range `r`.
"""
function diffraction_limited_spot_diameter(r::Unitful.Length, d_element::Unitful.Length, λ::Unitful.Length, m2 = 1)
	if r <= 0u"m"
		throw(DomainError(r, "Range r must be greater than zero"))
	elseif d_element <= 0u"m"
		throw(DomainError(d_element, "Element diameter must be greater than zero"))
	elseif λ <= 0u"m"
		throw(DomainError(λ, "Wavelength λ must be greater than zero"))
	elseif m2 < 1
		throw(DomainError(m2, "M² cannot be less than 1"))
	end

	if r < d_element
		@warn "Range r less than element size, result will be inaccurate" r element_size
	end

	return 2m2 * λ / (π * atan(d_element, 2r)) |> unit(d_element)
end

function diffraction_limited_element_size(r::Unitful.Length, d_spot::Unitful.Length, λ::Unitful.Length, m2 = 1)
	if r <= 0u"m"
		throw(DomainError(r, "Range r must be greater than zero"))
	elseif d_spot <= 0u"m"
		throw(DomainError(d_spot, "Spot diameter must be greater than zero"))
	elseif λ <= 0u"m"
		throw(DomainError(λ, "Wavelength λ must be greater than zero"))
	elseif m2 < 1
		throw(DomainError(m2, "M² cannot be less than 1"))
	end

	d_element =  4r * tan(m2 * λ / (π * d_spot))
	
	if r < d_element
		@warn "Range r less than element size, result will be inaccurate" r d_element
	end
	
	return d_element |> unit(d_spot)
end

function diffraction_limited_range(d_element::Unitful.Length, d_spot::Unitful.Length, λ::Unitful.Length, m2 = 1)
	if d_element <= 0u"m"
		throw(DomainError(d_element, "Element diameter must be greater than zero"))
	elseif d_spot <= 0u"m"
		throw(DomainError(d_spot, "Spot diameter must be greater than zero"))
	elseif λ <= 0u"m"
		throw(DomainError(λ, "Wavelength λ must be greater than zero"))
	elseif m2 < 1
		throw(DomainError(m2, "M² cannot be less than 1"))
	end

	r = d_element / 4tan(m2 * λ / (π * d_spot))
	
	if r < d_element
		@warn "Range r less than element size, result will be inaccurate" r d_element
	end
	
	return upreferred(r)
end

function diffraction_limited_wavelength(r::Unitful.Length, d_element::Unitful.Length, d_spot::Unitful.Length, m2 = 1)
	if d_element <= 0u"m"
		throw(DomainError(d_element, "Element diameter must be greater than zero"))
	elseif d_spot <= 0u"m"
		throw(DomainError(d_spot, "Spot diameter must be greater than zero"))
	elseif r <= 0u"m"
		throw(DomainError(r, "Range r must be greater than zero"))
	elseif m2 < 1
		throw(DomainError(m2, "M² cannot be less than 1"))
	end
	
	if r < d_element
		@warn "Range r less than element size, result will be inaccurate" r d_element
	end
	
	return tan(d_element /2r) * π * d_spot / 2m2 |> u"nm"
end

attenuated_beam_power(r::Unitful.Length, power::Unitful.Power, attenuation_length::Unitful.Length) = power * exp(-r / attenuation_length)

# http://panoptesv.com/SciFi/LaserDeathRay/Breakdown.html
function approximate_breakdown_threshold(λ::Unitful.Length, ambient_pressure::Unitful.Pressure)
	p = ustrip(ambient_pressure |> u"atm")
	l = ustrip(λ |> u"m")
	
	return ((0.31 / (l * p^0.6)^2) * u"W*cm^-2") |> u"W/m^2"
end

beam_intensity(power::Unitful.Power, d_spot::Unitful.Length) = power / (π * d_spot^2 / 4)

end