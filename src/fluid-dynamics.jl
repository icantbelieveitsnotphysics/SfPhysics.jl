module SfFluidDynamics

using Unitful, DifferentialEquations

import PhysicalConstants.CODATA2018: g_n

import ..SfPhysics: Angle, to_angle

export drag_force, reynolds_number

const Viscosity{T} = Unitful.AbstractQuantity{T, Unitful.𝐌 * Unitful.𝐋^-1 *  Unitful.𝐓^-1, typeof(u"Pa*s")}
	
"""
    drag_force(drag_coefficient::Real, ρ_fluid::Unitful.Density, reference_area::Unitful.Area, v::Unitful.Velocity)
	
Drag acting on a body with `drag_coefficient` and `reference_area` travelling at `v` through a fluid with density `ρ_fluid`.
"""
drag_force(drag_coefficient::Real, ρ_fluid::Unitful.Density, reference_area::Unitful.Area, v::Unitful.Velocity) = 
	upreferred(0.5ρ_fluid * v^2 * drag_coefficient * reference_area)

"""
    reynolds_number(v::Unitful.Velocity, ρ_fluid::Unitful.Density, l::Unitful.Length, μ::Viscosity)
	
Reynolds number of a body with characteristic dimension `l` travelling at velocity `v` through a fluid of density `ρ_fluid` and viscosity `μ`.
"""
reynolds_number(v::Unitful.Velocity, ρ_fluid::Unitful.Density, l::Unitful.Length, μ::Viscosity) = (v * ρ_fluid * l) / μ |>u"s/s"

# https://en.wikipedia.org/wiki/Projectile_motion#Numerical_solution

# 

"""
    ballistic_drag!(du, state, parameters, t)
    
simple projectile motion with drag ODE
 - const drag coefficient
 - constant air density
 - constant gravity

note that these things _could_ be implemented easily!
"""
function ballistic_drag!(du, state, parameters, t)
    # current state
    x, y, xv, yv = state
    # ballistic parameters
    μ, g = parameters
    # projectile velocity
    v = sqrt(xv^2 + yv^2)

    # d/dt
    du[1] = xv
    du[2] = yv
    du[3] = -μ * xv * v
    du[4] = -g - (μ * yv * v)
end

# note calls to ustrip... DifferentialEquations doesn't cope well with Unitful units, alas.
"""
    ballistic_initial_state(v0, θ0)
Initial state vector for `ballistic_drag!` simulation.

# Arguments
 - `v0`: Velocity magnitude. If not using units, this should be in metres per second.
 - `θ0`: Launch angle, assumed radians if no units given. 90° is straight up.
"""
ballistic_initial_state(v0, θ0) = [0.0;0.0;ustrip(cos(θ0)*(v0|>u"m/s"));ustrip(sin(θ0)*(v0|>u"m/s"))]

# parameters for simulation:
"""
    ballistic_params(Cd, ρ_a, A, m, g = g_n)

Parameters for `ballistic_drag!` simulation.

# Arguments
- `Cd`  : coefficient of drag for projectile
- `ρ_a` : density of atmosphere
- `A`   : cross-sectional area of projectile in direction of motion
- `m`   : projectile mass
- `g`   : acceleration due to gravity, defaults to Earth surface standard
"""
ballistic_params(Cd, ρ_a, A, m, g = g_n) = [ustrip((Cd*(ρ_a |> u"kg/m^3")*(A |> u"m^2")/2) / (m |> u"kg")), ustrip(g |> u"m/s^2")]

"""
    extract_velocity_vector_magnitude(u, vx, vy)

Velocity vector magnitude with arbitrary parameter as a 2-tuple, for plotting graphs using DifferentialEquations and Plots.
"""
extract_velocity_vector_magnitude(u, vx, vy) = (u, sqrt(vx^2 + vy^2))

end
