module SfPhysics

using Reexport

include("units.jl")

@reexport using .SfUnits

include("kinematics.jl")
include("fluid-dynamics.jl")

@reexport using .SfFluidDynamics

include("geometry.jl")

@reexport using .SfGeometry

include("gravity.jl")
include("blackhole.jl")

@reexport using .SfGravity
@reexport using .SfBlackHole

include("rocketry.jl")

@reexport using .SfRocketry

include("relativity.jl")

@reexport using .SfRelativity

include("matter.jl")

@reexport using .SfMatter

include("materials.jl")

@reexport using .SfMaterials

include("objects.jl")

@reexport using .SfObjects

include("astronomy.jl")

@reexport using .SfAstronomy

include("thermo.jl")

@reexport using .SfThermo

include("planetary.jl")
include("solar.jl")

@reexport using .SfPlanetary

export SfSolarSystem

include("orbital.jl")

@reexport using .SfOrbital

include("atmosphere.jl")

@reexport using .SfAtmosphere

include("laser.jl")

@reexport using .SfLaser

include("coriolis.jl")

@reexport using .SfCoriolis

include("electro.jl")

@reexport using .SfElectro

include("weaponry.jl")

@reexport using .SfWeaponry

end
