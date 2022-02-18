module SfAtmosphere

using Unitful
using ..SfPlanetary, ..SfThermo, ..SfMatter

import PhysicalConstants.CODATA2018: k_B, G, R, g_n
import ..SfUnits: MolarMass

export scale_height, jeans_escape_timescale, jeans_parameter

scale_height(t, m::Unitful.Mass, g = g_n) = (k_B * t) / (m * g) |> u"km"

scale_height(t, m::MolarMass, g = g_n) = (R * t) / (m * g) |> u"km"

function jeans_parameter(v_esc, v_peak)
    λ =(v_esc/v_peak)^2
    return ((1 + λ)*exp(-λ))/sqrt(4π)
end

function jeans_velocity(v_esc, m, T)
    v_peak=maxwell_boltzmann_peak_speed(m, T) # peak of maxwell-boltzmann distribution
    v_jeans = v_peak * jeans_parameter(v_esc, v_peak)
end

function maxwell_boltzmann_fraction(v, m, n = 2)
    
end

"""
    jeans_escape_timescale(T::Unitful.Temperature, exobase::Unitful.Length, M::Unitful.Mass, R::Unitful.Length, m::Unitful.Mass)
	
Approximate Jeans escape timescale for a gas of atomic mass `m` on a planet of mass `M` and radius `R`, with an exobase at altitude `xh` and temperature `T`.

Working from http://cococubed.asu.edu/code_pages/jeans_escape.shtm

"""
function jeans_escape_timescale(T::Unitful.Temperature, xh::Unitful.Length, M::Unitful.Mass, R::Unitful.Length, m::Unitful.Mass)
   g = (G*M)/(R+xh)^2 # gravity at exobase
   H = scale_height(T, m, g)
   v_esc=sqrt((2G*M)/(R+xh)) # escape velocity at exobase
   v_jeans = jeans_velocity(v_esc, m, T)
   return H/v_jeans |> u"yr"
end



"""
    jeans_escape_timescale(T::Unitful.Temperature, xh::Unitful.Length, body::AbstractBody, m::Unitful.Mass)
	
Approximate Jeans escape timescale for a gas of atomic mass `m` on planet `body`, with an exobase at altitude `xh` and temperature `T`.

Working from http://cococubed.asu.edu/code_pages/jeans_escape.shtm
"""
jeans_escape_timescale(T::Unitful.Temperature, xh::Unitful.Length, body::AbstractBody, m::Unitful.Mass) = 
	jeans_escape_timescale(T, xh, mass(body), radius(body), m)

# https://arxiv.org/ftp/arxiv/papers/1009/1009.5110.pdf
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008GL036513
jeans_parameter(m_planet::Unitful.Mass, m_molecule::Unitful.Mass, t_exosphere::Unitful.Temperature, r_exosphere::Unitful.Length) = 
	(G * m_planet * m_molecule) / (k_B * t_exosphere * r_exosphere) |> u"m/m"

jeans_parameter(body::AbstractBody, m_molecule::Unitful.Mass, t_exosphere::Unitful.Temperature, r_exosphere::Unitful.Length) = 
	jeans_parameter(mass(body), m_molecule, t_exosphere, r_exosphere)

end
