using PhysicalConstants.CODATA2014: c_0, g_n
using Unitful

function relativistic_brachistochrone_transit_time(dist::Unitful.Length, acc::Unitful.Acceleration)
    hd = dist/2

    2 * sqrt((hd/c_0)^2 + (2 * hd / acc))
end

function proper_relativistic_brachistochrone_transit_time(d::Unitful.Length, acc::Unitful.Acceleration)
    hd = d/2

    2 * (c_0 / acc) * acosh(((acc * hd) / c_0^2) + 1)
end

relativistic_velocity(t::Unitful.Time, acc::Unitful.Acceleration) = (acc*t)/sqrt(1+((acc*t)/c_0)^2|>u"s/s") |>u"m/s"

lorentz_factor(t::Unitful.Time, acc::Unitful.Acceleration) = sqrt(1 + (acc*t/c_0)^2) |> u"m/m"
lorentz_factor(d::Unitful.Length, acc::Unitful.Acceleration) = 1 + (acc * d)/c_0^2 |> u"m/m"

rocket_propulsive_efficiency(shipVelocity::Unitful.Velocity, exhaustVelocity::Unitful.Velocity) = 2(shipVelocity / exhaustVelocity) / (1+(shipVelocity / exhaustVelocity)^2)

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

relativistic_delta_v(ve::Unitful.Velocity, mr) = c_0 * tanh((ve/c_0) * log(mr)) |>u"m/s"

rocket_thrust(ve::Unitful.Velocity, mdot::Unitful.MassFlow) = Unitful.uconvert(Unitful.N, mdot * ve)
rocket_thrust(ve::Unitful.Velocity, fp::Unitful.Power) = Unitful.uconvert(Unitful.N, (2 * fp) / ve)

rocket_power(ve::Unitful.Velocity, mdot::Unitful.MassFlow) = Unitful.uconvert(Unitful.W, (mdot * ve^2) / 2)
rocket_power(ve::Unitful.Velocity, thrust::Unitful.Force) = Unitful.uconvert(Unitful.W, (thrust * ve) / 2)

mass_flow(ve::Unitful.Velocity, thrust::Unitful.Force) = Unitful.uconvert(u"kg/s", thrust / ve)

exhaust_velocity(isp::Unitful.Time) = isp * g_n
exhaust_velocity(thrust::Unitful.Force, mdot::Unitful.MassFlow) = Unitful.uconvert(u"m/s", thrust / mdot)
exhaust_velocity(fp::Unitful.Power, mdot::Unitful.MassFlow) = Unitful.uconvert(u"m/s", sqrt((2 * fp) / mdot))
