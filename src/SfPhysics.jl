module SfPhysics

using Reexport

include("kinematics.jl")
include("kinetics.jl")
include("fluid-dynamics.jl")

@reexport using .SfFluidDynamics

include("geometry.jl")
include("gravity.jl")
include("blackhole.jl")

@reexport using .SfGravity

include("rocketry.jl")

@reexport using .SfRocketry

include("relativity.jl")

@reexport using .SfRelativity

include("materials.jl")
include("planetary.jl")
include("solar.jl")

@reexport using .SfPlanetary

export SfSolarSystem

include("laser.jl")

@reexport using .SfLaser

include("coriolis.jl")

@reexport using .SfCoriolis

include("electro.jl")

@reexport using .SfElectro

include("weaponry.jl")

include("astronomy.jl")

@reexport using .SfAstronomy

end
