module SfVisualise

using Unitful

import ..SfOrbital: OrbitalElements, elements_to_state_vector, orbital_position, orbital_period

export xyz, quantise, sweep_quantise, plot_extent

xyz(v) = ([e[1] for e in v], [e[2] for e in v], [e[3] for e in v])

# evenly distributes samples around the orbit by true anomaly, but eccentric orbits will
# have more points near the focus and with a continuous time step they'll appear to move
# more slowly when you'd expect them to be moving fastest
function quantise(orbit::OrbitalElements, central_mass, steps)
    points = []
    
    o = deepcopy(orbit)
       
    for ν in 0:steps
       o.ν = deg2rad(ν)
       s = elements_to_state_vector(o, central_mass)
       push!(points, s.r)
    end
    
    return points
end

# even distribution of samples around the orbit by distance. this causes animations to
# look correct (fastest movement by the focus) but number of points near focus will be
# reduced for highly elliptic orbits which may impact fidelity
function sweep_quantise(orbit::OrbitalElements, central_mass, steps, time = :period; clip = false)
    points = []
    
    period = orbital_period(orbit, central_mass)
    
    if time == :period
        time = period
    end
    
    # the step period needs to be small for tight or highly elliptical orbits
    s = time / steps
    
    if clip && time > period
       steps =  period / s
    end
    
    t = 0u"s"
    
    for i in 0:steps
       o = orbital_position(orbit, t, central_mass)
       sv = elements_to_state_vector(o, central_mass)
       push!(points, sv.r)
       t += s
    end
    
    return points
end

function plot_extent(orbits)
    minval = typemin(orbits[1][1][1])
    maxval = typemax(orbits[1][1][1])
    
    xmin = ymin = zmin = maxval
    xmax = ymax = zmax = minval
    
    for orbit in orbits
        for point in orbit
            xmin = min(point[1], xmin)
            xmax = max(point[1], xmax)
            ymin = min(point[2], ymin)
            ymax = max(point[2], ymax)
            zmin = min(point[3], zmin)
            zmax = max(point[3], zmax)
        end
    end
    
    return ((xmin, xmax), (ymin, ymax), (zmin, zmax))
end

end