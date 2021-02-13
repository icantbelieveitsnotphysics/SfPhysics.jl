module SfPhysics

using Reexport

include("geometry.jl")
include("gravity.jl")
include("blackhole.jl")

@reexport using .SfGravity

include("rocketry.jl")
include("weaponry.jl")
include("kinematics.jl")
include("relativity.jl")

@reexport using .SfRelativity

include("planetary.jl")
include("solar.jl")

@reexport using .SfPlanetary

end
