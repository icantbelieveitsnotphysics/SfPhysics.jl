module SfRocketry

import PhysicalConstants.CODATA2018: c_0, g_n
using Unitful

import ..SfUnits: Angle, to_angle, Acceleration, Speed

export rocket_propulsive_efficiency, brachistochrone_transit_time, brachistochrone_acceleration, brachistochrone_delta_v,
	boost_coast_transit_time, boost_coast_thrust_time, beam_core_mass_ratio, mass_ratio, delta_v, rocket_thrust, rocket_power,
	mass_flow, exhaust_velocity, shkadov_thrust

rocket_propulsive_efficiency(ship_velocity::Speed, exhaust_velocity::Speed) = 
	2(ship_velocity / exhaust_velocity) / (1+(ship_velocity / exhaust_velocity)^2)

"""
    brachistochrone_transit_time(d::Unitful.Length, a::Acceleration)
	
Time to travel a distance of `d` given a constant acceleration of `a` and a final velocity of 0 relative to the start.
"""
brachistochrone_transit_time(d::Unitful.Length, a::Acceleration) = upreferred(2sqrt(d/a))

"""
	brachistochrone_acceleration(d::Unitful.Length, t::Unitful.Time)
	
Required constant acceleration to travel a distance of `d` in time `t` with a final velocity of 0 relative to the start.
"""
brachistochrone_acceleration(d::Unitful.Length, t::Unitful.Time) = upreferred(4d/t^2)

"""
	brachistochrone_delta_v(d::Unitful.Length, a::Acceleration)
	
Required ΔV to travel a distance `d` with constant acceleration `a` and attain a final velocity of 0 relative to the start.
"""
brachistochrone_delta_v(d::Unitful.Length, a::Acceleration) = 2sqrt(d*a) |> u"km/s"

"""
    boost_coast_transit_time(d::Unitful.Length, a::Acceleration, boost::Unitful.Time)
	
Total time to travel a distance of `d`, given speedup and slowdown phases of length `boost` undergoing acceleration `a`.
"""
boost_coast_transit_time(d::Unitful.Length, a::Acceleration, boost::Unitful.Time) = ((d - (a * boost^2)) / (a * boost)) + 2boost

"""
    boost_coast_thrust_time(dist::Unitful.Length, acc::Acceleration, transit::Unitful.Time)
	
Combined duration of boost and brake phases required to cross a distance of `dist` in time `transit` given an acceleration under thrust of `acc`.
"""
function boost_coast_thrust_time(dist::Unitful.Length, acc::Acceleration, transit::Unitful.Time)
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

"""
	beam_core_mass_ratio(a, dv::Speed, ve::Speed)
	
Mass ratio of a rocket where mass is lost without contributing to thrust, eg. a beam core rocket.

# Arguments
 - `a`: Fraction of propellant mass remaining after reaction (.22 for Frisbee rocket).
 - `dv::Speed`: mission delta-V (.25c for Frisbee rocket).
 - `ve::Speed`: exhaust velocity (.33c for Frisbee rocket).
 
 Derivation in Robert Frisbee's [How to build an antimatter rocket for interstellar missions](https://web.archive.org/web/20060601234257/http://www.aiaa.org/Participate/Uploads/2003-4676.pdf).
"""
function beam_core_mass_ratio(a, dv::Speed, ve::Speed)
	# un-tidied form of the equation looks more like this:
	# bcmr(a, ∆V, Isp) = ( ( (-2Isp*∆V/c_0^2+(1-a)-((1-a)^2+4a*Isp^2/c_0^2)^0.5)*(1-a+((1-a)^2+4a*Isp^2/c_0^2)^0.5) ) /( (-2Isp*∆V/c_0^2+(1-a)+((1-a)^2+4a*Isp^2/c_0^2)^0.5)*(1-a-((1-a)^2+4a*Isp^2/c_0^2)^0.5) ) )^(1 / ((1-a)^2+4a*Isp^2/c_0^2)^0.5)

	k1 = sqrt((1-a)^2+4a*ve^2/c_0^2)
	k2 = (-2ve * dv / c_0^2) + 1 - a

	p1 = k2 - k1
	p2 = 1 - a + k1
	p3 = k2 + k1
	p4 = 1 - a - k1
	p5 = 1 / k1

	((p1 * p2) / (p3 * p4))^p5
end

"""
    mass_ratio(dv::Speed, ve::Speed)
	
Compute the mass ratio for a rocket with total delta-V `dv` and exhaust velocity `ve`.
"""
mass_ratio(dv::Speed, ve::Speed) = exp(dv/ve)

delta_v(ve::Speed, mr::Real) = ve * log(mr)
delta_v(isp::Unitful.Time, mr::Real) = isp * g_n * log(mr)

rocket_thrust(ve::Speed, mdot::Unitful.MassFlow) = Unitful.uconvert(Unitful.N, mdot * ve)
rocket_thrust(ve::Speed, fp::Unitful.Power) = Unitful.uconvert(Unitful.N, (2 * fp) / ve)

rocket_power(ve::Speed, mdot::Unitful.MassFlow) = Unitful.uconvert(Unitful.W, (mdot * ve^2) / 2)
rocket_power(ve::Speed, thrust::Unitful.Force) = Unitful.uconvert(Unitful.W, (thrust * ve) / 2)

mass_flow(ve::Speed, thrust::Unitful.Force) = Unitful.uconvert(u"kg/s", thrust / ve)

"""
    exhaust_velocity(isp::Unitful.Time)
	
Convert specific impulse `isp` in seconds to exhaust velocity.
"""
exhaust_velocity(isp::Unitful.Time) = isp * g_n

"""
    exhaust_velocity(thrust::Unitful.Force, mdot::Unitful.MassFlow)
	
Compute exhaust velocity required for the given mass flow `mdot` to develop `thrust`.
"""
exhaust_velocity(thrust::Unitful.Force, mdot::Unitful.MassFlow) = Unitful.uconvert(u"m/s", thrust / mdot)

"""
    exhaust_velocity(fp::Unitful.Power, mdot::Unitful.MassFlow)
	
Compute exhaust velocity required for the given mass flow `mdot` and thrust power `fp`.
"""
exhaust_velocity(fp::Unitful.Power, mdot::Unitful.MassFlow) = Unitful.uconvert(u"m/s", sqrt((2 * fp) / mdot))

"""
    shkadov_thrust(stellar_luminosity::Unitful.Power, rim_angle::Angle)
	
Compute thrust of a Shkadov thruster with given rim angle (90° for max thrust_ about a star with the given `stellar_luminosity` in watts.

Note that this is a simple photon thruster.
"""
shkadov_thrust(stellar_luminosity::Unitful.Power, rim_angle::Angle) = (stellar_luminosity / 2c_0) * (1 - cos(rim_angle)) |> u"N"

end
