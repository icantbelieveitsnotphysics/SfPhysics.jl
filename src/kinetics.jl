using Unitful

export kinetic_energy

"""
	kinetic_energy(m::Unitful.Mass, v::Unitful.Velocity)
	
Compute the kinetic energy of a body with mass `m` travelling at velocity `v`.

No relativistic corrections are applied. Use `relativistic_kinetic_energy` if those are required.
"""
kinetic_energy(m::Unitful.Mass, v::Unitful.Velocity) = 0.5m*v^2 |> u"J"