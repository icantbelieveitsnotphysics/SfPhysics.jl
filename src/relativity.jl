module SfRelativity

	using PhysicalConstants.CODATA2018: c_0, g_n, G
	using Unitful
	using UnitfulAstro

	#import ..
	
	export lorentz_factor, lorentz_velocity, relativistic_kinetic_energy

	lorentz_factor(v::Unitful.Velocity) = 1/sqrt(1-(v/c_0)^2)
	lorentz_velocity(γ) = sqrt(1 - (1/γ)^2) * c_0 |> u"c"

	function relativistic_kinetic_energy(m::Unitful.Mass, v::Unitful.Velocity)
		if v < 100u"km/s"
			return kinetic_energy(m, v)
		end

		e = m * big(c_0)^2

		e * (lorentz_factor(v) - 1) |> u"J"
	end
	
end