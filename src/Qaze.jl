module Qaze
using TOML
include("Structures.jl")
include("Config.jl")
include("BlackHole.jl")
include("Constants.jl")
include("InitialConditions.jl")
include("Integrate.jl")
#include("Plotting.jl")
include("QuadTree.jl")
include("QuadTreeInter.jl")
include("Radiation.jl")
include("Streamline.jl")
include("Utils.jl")
include("Wind.jl")

end #module
