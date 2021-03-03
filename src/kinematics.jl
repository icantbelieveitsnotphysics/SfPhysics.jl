using Unitful

export distance, duration, acceleration

"""
    distance(a::Unitful.Acceleration, t::Unitful.Time, initial_v::Unitful.Velocity = 0u"m/s")
	
Compute the distance travelled in time `t` by a body with initial velocity `initial_v` and uniform acceleration `a`.
"""
distance(a::Unitful.Acceleration, t::Unitful.Time, initial_v::Unitful.Velocity = 0u"m/s") = t*initial_v + a*t^2 / 2

"""
    duration(a::Unitful.Acceleration, d::Unitful.Length, initial_v::Unitful.Velocity = 0u"m/s")
	
Compute the time taken for a body with initial velocity `initial_v` and uniform acceleration `a` to travel a distance of `d`.
"""
function duration(a::Unitful.Acceleration, d::Unitful.Length, initial_v::Unitful.Velocity = 0u"m/s")
    if d < 0u"m"
		throw(DomainError(d, "Negative distances not allowed."))
	end

	if initial_v == 0u"m/s"
		return sqrt(2d / a)
	elseif a == 0u"m/s^2"
		if initial_v < 0u"m/s"
			throw(DomainError(d, "Negative initial velocities not allowed with zero acceleration."))
		end
	
		return d / initial_v
	end

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

"""
    acceleration(d::Unitful.Length, v_final::Unitful.Velocity, v_initial::Unitful.Velocity = 0u"m/s")
	
Given a distance `d`, initial velocity `v_initial and desired velocity `v_final` compute the required uniform acceleration
"""
acceleration(d::Unitful.Length, v_final::Unitful.Velocity, v_initial::Unitful.Velocity = 0u"m/s") = (v_final^2 - v_initial^2) / 2d