module SfAstronomy

using Unitful, UnitfulAstro

import ..SfUnits: Angle, to_angle

export absolute_magnitude, apparent_magnitude, diffuse_disc_q, diffuse_sphere_q

"""
    absolute_magnitude(radius::Unitful.Length, geometric_albedo::Real)
	
Compute the absolute magnitude of lambertian disc with the given `radius` and `geometric_albedo`.
"""
function absolute_magnitude(radius::Unitful.Length, geometric_albedo::Real)
	d = ustrip(radius *2 |> u"km")
	return 5 * log10(1329 / (d * sqrt(geometric_albedo)))
end

"""
	diffuse_sphere_q(α::Angle)
	
Approximate phase integral of an ideal diffuse reflecting sphere for observer angle `α`.
"""
diffuse_sphere_q(α::Angle) = (2/3) * ((1 - (α / 180u"°")) * cos(α) + (1/π) * sin(α))

"""
	diffuse_disc_q(α::Angle)
	
Phase integral of a diffuse reflecting disc for observer angle `α`.
"""
diffuse_disc_q(α::Angle) = cos(α)

"""
	apparent_magnitude(h::Real, d_bs::Unitful.Length, d_bo::Unitful.Length, d_os::Unitful.Length, q_α::Real)
	
Compute the apparent magnitude of an object with absolute magnitude `h`, distance to the sun `d_bs`, distance to the observer `d_bo`, distance from the observer to the sun `d_os` and phase integral `q_α`.
"""
function apparent_magnitude(h::Real, d_bs::Unitful.Length, d_bo::Unitful.Length, d_os::Unitful.Length, q_α::Real)
	return h + 5 * log10(ustrip(d_bs * d_bo / d_os^2 |> u"km/km")) - 2.5 * log10(q_α)
end

end
