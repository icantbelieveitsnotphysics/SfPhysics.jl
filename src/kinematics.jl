﻿using Unitful

import PhysicalConstants.CODATA2018: g_n

export distance, duration, acceleration, projectile_displacement, projectile_velocity, projectile_flight_time, projectile_peak_displacement

# making vectors of abstract unitful types like Unitful.Acceleration is a right pain, and this is much simpler
const Accel{T} = Unitful.AbstractQuantity{T,Unitful.𝐋*Unitful.𝐓^-2,typeof(u"m/s/s")}
const AccelerationVector{T} = Vector{Accel{T}}
const Speed{T} = Quantity{T,Unitful.𝐋*Unitful.𝐓^-1,typeof(u"m/s")}
const SpeedVector{T} = Vector{Speed{T}}

const Angle{T} = Union{ Quantity{T, NoDims, typeof(u"°")}, Quantity{T, NoDims, typeof(u"rad")} }

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

"""
    projectile_displacement(v_0::Unitful.Velocity, θ::Angle, t::Unitful.Time, g::Accel = g_n)
	
Horizontal and vertical displacement at time `t` of a projectile with initial velocity `v_0` and launch angle `θ` in a unitform gravitational acceleration `g`.
"""
projectile_displacement(v_0::Unitful.Velocity, θ::Angle, t::Unitful.Time, g::Accel = g_n) =
    [ v_0 * t * cos(θ), v_0 * t * sin(θ) - 0.5g * t^2 ]

"""
    projectile_displacement(v_0::Unitful.Velocity, t::Unitful.Time, g::Accel = g_n)
	
Vertical displacement at time `t` of a projectile launched vertically with initial velocity `v_0` in a unitform gravitational acceleration `g`.
"""
projectile_displacement(v_0::Unitful.Velocity, t::Unitful.Time, g::Accel = g_n) = projectile_displacement(v_0, 90u"°", g, t)[2]

"""
    projectile_velocity(v_0::Unitful.Velocity, θ::Angle, t::Unitful.Time, g::Accel = g_n)
	
Horizontal and vertical velocity at time `t` of a projectile with initial velocity `v_0` and launch angle `θ` in a unitform gravitational acceleration `g`.
"""
projectile_velocity(v_0::Unitful.Velocity, θ::Angle, t::Unitful.Time, g::Accel = g_n) =
    [ v_0 * cos(θ), v_0 * sin(θ) - g * t ]

"""
    projectile_velocity(v_0::Unitful.Velocity, t::Unitful.Time, g::Accel = g_n)
	
Velocity at time `t` of a projectile launched vertically with initial velocity `v_0` in a unitform gravitational acceleration `g`.
"""
projectile_velocity(v_0::Unitful.Velocity, t::Unitful.Time, g::Accel = g_n) = projectile_velocity(v_0, 90u"°", g, t)[2]

"""
    projectile_flight_time(v_0::Unitful.Velocity, θ::Angle = 90u"°", g::Accel = g_n)
	
Compute total flight time for a projectile launched at angle `θ` with initial velocity `v_0` in a uniform gravitation acceleration `g`.

Flight is considered to be complete when vertical displacement returns to zero.
"""
projectile_flight_time(v_0::Unitful.Velocity, θ::Angle = 90u"°", g::Accel = g_n) = 2v_0 * sin(θ) / g |> u"s"

"""
    projectile_peak_displacement(v_0::Unitful.Velocity, θ::Angle = 90u"°", g::Accel = g_n)
	
Compute peak altitude of a projectile launched at angle `θ` with initial velocity `v_0` in a uniform gravitation acceleration `g`.
"""
projectile_peak_displacement(v_0::Unitful.Velocity, θ::Angle = 90u"°", g::Accel = g_n) = v_0^2 / 2g