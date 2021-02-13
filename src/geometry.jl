using Unitful

sphere_volume(r::Unitful.Length)::Unitful.Volume = (4π/3)r^3
sphere_radius(v::Unitful.Volume)::Unitful.Length = cbrt(3v / 4π)

cylinder_volume(r::Unitful.Length, h::Unitful.Length)::Unitful.Volume = π * r^2 * h
cylinder_radius(v::Unitful.Volume, h::Unitful.Length)::Unitful.Length = sqrt(v / (π * h))
cylinder_length(v::Unitful.Volume, r::Unitful.Length)::Unitful.Length = v / (π * r^2)
