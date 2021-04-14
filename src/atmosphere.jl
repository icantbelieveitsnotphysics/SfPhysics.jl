module SfAtmosphere

using Unitful
using ..SfPlanetary, ..SfThermo, ..SfMatter

import PhysicalConstants.CODATA2018: k_B, G, R, g_n
import ..SfUnits: MolarMass

export scale_height, jeans_escape_timescale, jeans_parameter

scale_height(t, m::Unitful.Mass, g = g_n) = (k_B * t) / (m * g) |> u"km"

scale_height(t, m::MolarMass, g = g_n) = (R * t) / (m * g) |> u"km"

# temperature, exosphere altiutude, planetary mass, planetary radius, gas molecular mass
# http://cococubed.asu.edu/code_pages/jeans_escape.shtml
function jeans_escape_timescale(T::Unitful.Temperature, h::Unitful.Length, M::Unitful.Mass, R::Unitful.Length, m::Unitful.Mass)
   g = (G*M)/(R+h)^2 # gravity at exobase
   H = (k_B*T)/(m*g) # scale height for gas
   v_peak=maxwell_boltzmann_peak_speed(m, T) # peak of maxwell-boltzmann distribution
   v_esc=sqrt((2G*M)/(R+h)) # escape velocity of planet at exosphere altitude
   λ =(v_esc/v_peak)^2
   v_jeans = v_peak * ((1 + λ)*exp(-λ))/sqrt(4π)
   return H/v_jeans |> u"yr"
end

jeans_escape_timescale(T::Unitful.Temperature, h::Unitful.Length, body::AbstractBody, m::Unitful.Mass) = 
	jeans_escape_timescale(T, h, mass(body), radius(body), m)

# https://arxiv.org/ftp/arxiv/papers/1009/1009.5110.pdf
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008GL036513
jeans_parameter(m_planet::Unitful.Mass, m_molecule::Unitful.Mass, t_exosphere::Unitful.Temperature, r_exosphere::Unitful.Length) = 
	(G * m_planet * m_molecule) / (k_B * t_exosphere * r_exosphere) |> u"m/m"

jeans_parameter(body::AbstractBody, m_molecule::Unitful.Mass, t_exosphere::Unitful.Temperature, r_exosphere::Unitful.Length) = 
	jeans_parameter(mass(body), m_molecule, t_exosphere, r_exosphere)

end
