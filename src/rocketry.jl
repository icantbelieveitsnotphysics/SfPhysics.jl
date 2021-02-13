module SfRocketry

	using PhysicalConstants.CODATA2018: c_0, g_n
	using Unitful
	
	export rocket_propulsive_efficiency, brachistochrone_transit_time, brachistochrone_acceleration, brachistochrone_delta_v,
		boost_coast_transit_time, coasting_coast_time, beam, mass_ratio, delta_v, rocket_thrust, rocket_power,
		mass_flow, exhaust_velocity

	rocket_propulsive_efficiency(ship_velocity::Unitful.Velocity, exhaust_velocity::Unitful.Velocity) = 2(ship_velocity / exhaust_velocity) / (1+(ship_velocity / exhaust_velocity)^2)

	brachistochrone_transit_time(d::Unitful.Length, a::Unitful.Acceleration) = 2sqrt(d/a)
	brachistochrone_acceleration(d::Unitful.Length, t::Unitful.Time) = 4d/t^2
	brachistochrone_delta_v(d::Unitful.Length, a::Unitful.Acceleration) = 2sqrt(d*a)

	boost_coast_transit_time(d::Unitful.Length, a::Unitful.Acceleration, boost::Unitful.Time) = ((d - (a * boost^2)) / (a * boost)) + 2boost
	function coasting_coast_time(dist::Unitful.Length, acc::Unitful.Acceleration, transit::Unitful.Time)
		# Tt - t^2 - d/a
		# ax^2 + bx + c

		a = -1
		b = transit
		c = -dist / acc

		discr = b^2 - 4*a*c

		t1 = (-b + sqrt(discr))/(2a)
		t2 = (-b - sqrt(discr))/(2a)

		if t1 < 0u"s"
			return t2
		elseif t2 < 0u"s"
			return t1
		else
			return min(t1, t2)
		end
	end

	function beam_core_mass_ratio(a, dv::Unitful.Velocity, ve::Unitful.Velocity)
		k1 = sqrt((1-a)^2+4a*ve^2/c_0^2)
		k2 = (-2ve * dv / c_0^2) + 1 - a

		p1 = k2 - k1
		p2 = 1 - a + k1
		p3 = k2 + k1
		p4 = 1 - a - k1
		p5 = 1 / k1

		((p1 * p2) / (p3 * p4))^p5
	end

	mass_ratio(dv::Unitful.Velocity, ve::Unitful.Velocity) = exp(dv/ve)

	delta_v(ve::Unitful.Velocity, mr) = ve * log(mr)
	delta_v(isp::Unitful.Time, mr) = isp * g_n * log(mr)

	rocket_thrust(ve::Unitful.Velocity, mdot::Unitful.MassFlow) = Unitful.uconvert(Unitful.N, mdot * ve)
	rocket_thrust(ve::Unitful.Velocity, fp::Unitful.Power) = Unitful.uconvert(Unitful.N, (2 * fp) / ve)

	rocket_power(ve::Unitful.Velocity, mdot::Unitful.MassFlow) = Unitful.uconvert(Unitful.W, (mdot * ve^2) / 2)
	rocket_power(ve::Unitful.Velocity, thrust::Unitful.Force) = Unitful.uconvert(Unitful.W, (thrust * ve) / 2)

	mass_flow(ve::Unitful.Velocity, thrust::Unitful.Force) = Unitful.uconvert(u"kg/s", thrust / ve)

	exhaust_velocity(isp::Unitful.Time) = isp * g_n
	exhaust_velocity(thrust::Unitful.Force, mdot::Unitful.MassFlow) = Unitful.uconvert(u"m/s", thrust / mdot)
	exhaust_velocity(fp::Unitful.Power, mdot::Unitful.MassFlow) = Unitful.uconvert(u"m/s", sqrt((2 * fp) / mdot))

end