module SfOrbital

using Unitful, UnitfulAstro, UnitfulAngles
import PhysicalConstants.CODATA2018: G
import LinearAlgebra: norm, dot, cross

using ..SfGeometry
using ..SfPlanetary

import ..SfUnits: Angle, to_angle
import ..SfGravity: orbital_radius, orbital_period

export OrbitalElements, StateVector, state_vector_to_elements, elements_to_state_vector, time_since_periapsis,
    orbital_position

# very useful resource: https://orbital-mechanics.space

mutable struct OrbitalElements
    e # eccentricity
    a # semi-major axis
    i # inclination
    Ω # longitude of ascending node
    ω # argument of periapsis
    ν # true anomaly
end

OrbitalElements(;e=0, a, i=0, Ω=0, ω=0, ν=0) = OrbitalElements(e, a, i, Ω, ω, ν)

struct StateVector
    r # position
    v # velocity
end

function state_vector_to_elements(orbit::StateVector, central_mass)
    r = norm(orbit.r)
    v = norm(orbit.v)
    v_r = dot(r/norm(r), v) # radial velocity
    v_⊥ = sqrt(v^2 - v_r^2) # azimuthal velocity
    
    h_vec = cross(orbit.r, orbit.v)
    h = norm(h_vec) # orbital angular momentum
    
    i = acos(h_vec[3] / h) # inclination
    
    N_vec = cross([0, 0, 1], h_vec)
    N = norm(N_vec)
    Ω = 2π - acos(N_vec[1] / N) # right ascension
    
    μ = G * central_mass
    e_vec = cross(orbit.v, h_vec) / μ - orbit.r / r
    e = norm(e_vec) # eccentricity
    
    ω = 2π - acos(dot(N_vec, e_vec) / (N * e)) # argument of periapsis

    if v_r < 0u"m/s"
        ν = 2π - acos(dot(e_vec, orbit.r) / (e * r)) # true anomaly
    else
        ν = acos(dot(e_vec, orbit.r) / (e * r)) # true anomaly
    end
    
    # https://en.wikipedia.org/wiki/Semi-major_and_semi-minor_axes#Orbital_period
    a = h^2 / (μ * (1 - e^2)) # semi-major axis
    
    OrbitalElements(e, upreferred(a), i * u"rad", Ω * u"rad", ω * u"rad", (ν % 2π) * u"rad")
end

function elements_to_state_vector(orbit::OrbitalElements, central_mass)
    μ = G * central_mass
    
    h = sqrt(orbit.a * μ * (1 - orbit.e^2)) # orbital angular momentum
    
    # transform to perifocal plane
    r_w = h^2 / μ / (1 + orbit.e * cos(orbit.ν)) .* [cos(orbit.ν) sin(orbit.ν) 0]
    v_w = μ / h .* [-sin(orbit.ν) orbit.e + cos(orbit.ν) 0]
    
    # rotate perifocal plane to inertial plane
    R1 = [cos(-orbit.ω) -sin(-orbit.ω) 0; sin(-orbit.ω) cos(-orbit.ω) 0; 0 0 1]
    R2 = [1 0 0; 0 cos(-orbit.i) -sin(-orbit.i); 0 sin(-orbit.i) cos(-orbit.i)];
    R3 = [cos(-orbit.Ω) -sin(-orbit.Ω) 0; sin(-orbit.Ω) cos(-orbit.Ω) 0; 0 0 1]
    r_rot = r_w * R1 * R2 * R3
    v_rot = v_w * R1 * R2 * R3
    
    StateVector(upreferred.(vec(r_rot)), upreferred.(vec(v_rot)))
end

mean_motion(orbit::OrbitalElements, central_mass) = (2π * u"rad") / orbital_period(orbit, central_mass)

"""
    orbital_position(orbit::OrbitalElements, t, central_mass, ϵ = 0.000001)
    
Given some Keplerian elements for an object of negligible mass orbiting a body with mass `central_mass`,
compute its orbital elements after some elapsed time `t`.
"""
function orbital_position(orbit::OrbitalElements, t, central_mass, ϵ = 0.000001)
    # mutable struct, so always deep copy. 
    # all elements but true anomaly will be the same in the result
    pos = deepcopy(orbit)
    
    if orbit.ν != 0u"rad"
        t += time_since_periapsis(orbit, central_mass)
    end
    
    M = t * mean_motion(orbit, central_mass) # mean anomaly
    
    # circular orbits are simple
    if orbit.e == 1
        pos.ν = M
        return pos
    end
    
    if ustrip(t) == 0 || ustrip(t % orbital_period(orbit, central_mass)) < ϵ # close enough to periapsis
        pos.ν = 0u"rad"
        return pos
    end
    
    # newton's method
    f(est) = est - orbit.e * sin(est) - M
    f′(est) = 1 - orbit.e * cos(est)
    
    # sensible initial value
    if orbit.e > 0.8
        E = M
    else
        E = π
    end
    
    # iterate until a good enough approximation of the eccentric anomaly is found
    while abs(f(E)) > ϵ
        E = E - (f(E) / f′(E))
    end
    
    # true anomaly
    # this method is safe from quadrant problems
    # https://orbital-mechanics.space/time-since-periapsis-and-keplers-equation/elliptical-orbits.html
    pos.ν = 2atan(tan(E / 2) / sqrt((1 - orbit.e)/(1 + orbit.e)))
    
    return pos
end

"""
    time_since_periapsis(orbit::OrbitalElements, central_mass)
    
Compute the time elapsed for a massless body orbiting a body of mass `central_mass`
to have moved from its periapsis to its current true anomaly.
"""
function time_since_periapsis(orbit::OrbitalElements, central_mass)
    # circular orbits are simple
    if orbit.e == 1
        M = orbit.ν
    else
        c = cos(orbit.ν)
        E = acos((orbit.e + c)/(1 + orbit.e * c)) # eccentric anomaly
        M = E - orbit.e * sin(E) # mean anomaly
    end
    
    n = mean_motion(orbit, central_mass)
    
    upreferred(M / n)
end

"""
	orbital_radius(orbit::OrbitalElements)
	
Calculate the current distance of a body with Keplerian elements `orbit` from the body it is orbiting.
"""
orbital_radius(orbit::OrbitalElements) = orbit.a * ((1-orbit.e^2) / (1 + orbit.e * cos(orbit.ν)))

"""
    orbital_period(orbit::OrbitalElements, central_mass)

Calculate the orbital period of a body with Keplerian elements `orbit`, orbiting
a body with mass `central_mass`.
"""
orbital_period(orbit::OrbitalElements, central_mass) = orbital_period(central_mass, orbit.a)

end