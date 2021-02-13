module SfGravity

	using PhysicalConstants.CODATA2014: c_0, g_n, G, StefanBoltzmannConstant, ħ, k_B
	using Unitful
	using UnitfulAstro
	
	export vis_viva, gravity, planetary_mass, planetary_radius, escape_velocity, hill_sphere

	function vis_viva(parent_mass::Unitful.Mass, semimajor_axis::Unitful.Length, current_radius::Unitful.Length)
		G*parent_mass*(2/current_radius - 1/semimajor_axis) |> u"m^2/s^2"
	end

	gravity(m::Unitful.Mass, r::Unitful.Length) = G * m / r^2 |> u"m/s/s"
	gravity(m1::Unitful.Mass, m2::Unitful.Mass, r::Unitful.Length) = G * m * m / r^2 |> u"N"
	planetary_mass(acc::Unitful.Acceleration, r::Unitful.Length) = (acc * r^2) / G |> u"kg"
	planetary_radius(m::Unitful.Mass, acc::Unitful.Acceleration) = sqrt((G * m) / acc) |> u"m"

	# kepler's third
	orbital_period(m::Unitful.Mass, r::Unitful.Length) = sqrt((4π^2 * r^3) / (G * m)) |> u"s"
	orbital_radius(m::Unitful.Mass, t::Unitful.Time) = cbrt((t^2 * G * m)/(4π^2)) |> u"m"
	planetary_mass(orbital_radius::Unitful.Length, orbital_period::Unitful.Time) = (4π^2 * orbital_radius^3) / (G * orbital_period^2) |> u"kg"

	escape_velocity(m::Unitful.Mass, r::Unitful.Length) = sqrt(2G * m / r) |> u"km/s"
	planetary_mass(r::Unitful.Length, v_esc::Unitful.Velocity) = v_esc^2 * r / 2G |> u"kg"
	planetary_radius(m::Unitful.Mass, v_esc::Unitful.Velocity) = 2G * m / v_esc^2 |> u"m"

	function hill_sphere(m_parent::Unitful.Mass, m::Unitful.Mass, sma::Unitful.Length, e)
		sma * (1-e) * cbrt(m / 3m_parent)
	end
	hill_sphere(m_parent::Unitful.Mass, m::Unitful.Mass, sma::Unitful.Length) = sma * cbrt(m / 3m_parent) |> Unitful.unit(sma)

end