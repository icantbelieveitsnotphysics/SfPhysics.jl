using Unitful

export drag_force
	
"""
    drag_force(drag_coefficient::Real, ρ_fluid::Unitful.Density, reference_area::Unitful.Area, v::Unitful.Velocity)
	
Drag acting on a body with `drag_coefficient` and `reference_area` travelling at `v` through a fluid with density `ρ_fluid`.
"""
drag_force(drag_coefficient::Real, ρ_fluid::Unitful.Density, reference_area::Unitful.Area, v::Unitful.Velocity) = 
	upreferred(0.5ρ_fluid * v^2 * drag_coefficient * reference_area)