using Unitful, UnitfulAstro

import PhysicalConstants.CODATA2018: g_n

import ..SfUnits: Angle, to_angle

export kinetic_energy, distance, duration, acceleration, projectile_displacement, projectile_velocity, 
	projectile_flight_time, projectile_peak_displacement, projectile_range, projectile_angle, projectile_angle_planetary
	
# TODO: rename these when the underlying compilation issues have been resolved
export projectile_range_planetary
# TODO: better names?
export projectile_peak_displacement_planetary

"""
	kinetic_energy(m::Unitful.Mass, v::Unitful.Velocity)
	
Compute the kinetic energy of a body with mass `m` travelling at velocity `v`.

No relativistic corrections are applied. Use `relativistic_kinetic_energy` if those are required.
"""
kinetic_energy(m::Unitful.Mass, v::Unitful.Velocity) = 0.5m*v^2 |> u"J"

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

Default launch angle is 90° for maximum height.
"""
projectile_peak_displacement(v_0::Unitful.Velocity, θ::Angle = 90u"°", g::Accel = g_n) = (v_0^2 * sin(θ)^2) / 2g

"""
    projectile_range(v_0::Unitful.Velocity, θ::Angle = 45u"°", g::Accel = g_n)
	
Compute range of a projectile launched at angle `θ` with initial velocity `v_0` in a uniform gravitation acceleration `g`.

Default launch angle is 45° for maximum range.
"""
projectile_range(v_0::Unitful.Velocity, θ::Angle = 45u"°", g::Accel = g_n) = (v_0^2 * sin(2θ)) / g

"""
    projectile_range(v_0::Unitful.Velocity, y_0::Unitful.Length, θ::Angle = 45u"°", g::Accel = g_n)
	
Compute range of a projectile launched from altitude `y_0` at angle `θ` with initial velocity `v_0` in a uniform gravitation acceleration `g`.

Range is reached when vertical displacement is zero. Default launch angle is 45° for maximum range.
"""
projectile_range(v_0::Unitful.Velocity, y_0::Unitful.Length, θ::Angle = 45u"°", g::Accel = g_n) =
	(v_0 * cos(θ) / g) * (v_0 * sin(θ) + sqrt((v_0 * sin(θ))^2 + 2g * y_0))
	
"""
    projectile_range_planetary(v_0::Unitful.Velocity, r_planet::Unitful.Length = 1u"Rearth", θ::Angle = 45u"°", g::Accel = g_n)
	
Compute range of a projectile launched from altitude `y_0` at angle `θ` with initial velocity `v_0` on a spherical planet with radius `r_planet` and surface gravity `g`.

NOTE: this should be called `projectile_range`, but weird internal compiler errors prevent that with this version of Julia (1.5.3).
"""
function projectile_range_planetary(v_0::Unitful.Velocity, θ::Angle = 45u"°", r_planet::Unitful.Length = 1u"Rearth", g::Accel = g_n)
	v_rat2 = (v_0 / sqrt(r_planet * g))^2
	
	if v_rat2 > 1
		throw(DomainError("Initial velocity exceeds orbital velocity; Range undefined"))
	end
	
	a = v_0^2 * sin(2θ) / g
	b = sqrt(1 - (2 - v_rat2) * v_rat2 * cos(θ)^2)
	
	return a / b |> u"m"
end

function projectile_peak_displacement_planetary(v_0::Unitful.Velocity, θ::Angle = 45u"°", r_planet::Unitful.Length = 1u"Rearth", g::Accel = g_n)
	v_rat2 = (v_0 / sqrt(r_planet * g))^2
	
	if v_rat2 > 1
		throw(DomainError("Initial velocity exceeds orbital velocity; Range undefined"))
	end
	
	a = v_0^2 * sin(θ) / g
	b = 1 - v_rat2 + sqrt(1 - (2 - v_rat2) * v_rat2 * cos(θ)^2)
	
	return a / b |> u"m"
end
	
"""
    projectile_angle(v_0::Unitful.Velocity, d::Unitful.Length, g::Accel = g_n)
	
Compute the possible launch angles for a projectile to reach a horizontal range of `d` given initial velocity `v_0` in a uniform gravitation acceleration `g`.

The shallow angle is the first result, the steep angle the second. A domain error is raised if the target is out of range.
"""
function projectile_angle(v_0::Unitful.Velocity, d::Unitful.Length, g::Accel = g_n)
	k = (g * d) / v_0^2
	
	if (k > 1)
		throw(DomainError("Target out of range; no solutions"))
	end
	
	return ( 0.5asin(k) * 1u"rad" |> u"°", 45u"°" + (0.5acos(k) * 1u"rad" |> u"°") )
end
	
"""
    projectile_angle(v_0::Unitful.Velocity, d::Unitful.Length, y::Unitful.Length, g::Accel = g_n)
	
Compute the possible launch angles for a projectile to reach a horizontal range of `d` and altitude of `y` relative to the starting point given initial velocity `v_0` in a uniform gravitation acceleration `g`.

The shallow angle is the first result, the steep angle the second. A domain error is raised if the target is out of range.
"""
function projectile_angle(v_0::Unitful.Velocity, d::Unitful.Length, y::Unitful.Length, g::Accel = g_n)
	det = v_0^4 - g * (g * d^2 + 2y * v_0^2)
	
	if (ustrip(det) < 0)
		throw(DomainError("Target out of range; no solutions"))
	end
	
	a = atan((v_0^2 + sqrt(det)) / (g * d)) * 1u"rad" |> u"°"
	b = atan((v_0^2 - sqrt(det)) / (g * d)) * 1u"rad" |> u"°"
	
	if a < b
		return (a, b)
	else
		return (b, a)
	end
end

"""
    projectile_optimum_angle(v_0::Unitful.Velocity, r_planet::Unitful.Length = 1u"Rearth", g::Accel = g_n)
"""
function projectile_optimum_angle(v_0::Unitful.Velocity, r_planet::Unitful.Length = 1u"Rearth", g::Accel = g_n)
	v_rat2 = (v_0 / sqrt(r_planet * g))^2
	
	if v_rat2 > 0
		throw(DomainError("Initial velocity exceeds orbital velocity; optimum angle undefined"))
	end
	
	return 0.5acos(v_rat2 / (2 - v_rat2)) * u"rad" |> u"°"
end
