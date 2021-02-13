using Unitful

export distance, duration, kinetic_energy

distance(a::Unitful.Acceleration, t::Unitful.Time, initial_v::Unitful.Velocity) = t*initial_v + a*t^2 / 2
distance(a::Unitful.Acceleration, t::Unitful.Time)::Unitful.Length = distance(a, t, 0u"m/s")

duration(a::Unitful.Acceleration, d::Unitful.Length)::Unitful.Time = sqrt(2d / a)
function duration(a::Unitful.Acceleration, d::Unitful.Length, initial_v::Unitful.Velocity)::Unitful.Time
    # at^2 + 2v_0t - 2d = 0
    # ax^2 + bx + c

    b = 2 * initial_v
    c = -2d

    discr = b^2 - 4*a*c

    t1 = (-b + sqrt(discr))/(2a)
    t2 = (-b - sqrt(discr))/(2a)

    # the negative results aren't wrong, per se, but as the parameters are
    # vector quantities the negative results aren't particularly interesting.    

    if t1 < 0u"s"
        return t2
    elseif t2 < 0u"s"
        return t1
    else
        return min(t1, t2)
    end
end

kinetic_energy(m::Unitful.Mass, v::Unitful.Velocity) = 0.5m*v^2 |> u"J"
