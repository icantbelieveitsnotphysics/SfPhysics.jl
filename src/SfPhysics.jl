module SfPhysics

using Reexport

include("geometry.jl")
include("gravity.jl")
include("blackhole.jl")

@reexport using .SfGravity

include("rocketry.jl")

@reexport using .SfRocketry
include("kinematics.jl")
include("kinetics.jl")
include("relativity.jl")

@reexport using .SfRelativity

include("planetary.jl")
include("solar.jl")

@reexport using .SfPlanetary

export SfSolarSystem

include("materials.jl")
include("laser.jl")

@reexport using .SfLaser

include("coriolis.jl")

@reexport using .SfCoriolis

include("electro.jl")

@reexport using .SfElectro

include("weaponry.jl")

end
