module Qaze
export initialize_line, start_lines
using TOML
include("Structures.jl")
include("Config.jl")
include("BlackHole.jl")
include("Constants.jl")
include("InitialConditions.jl")
include("Integrate.jl")
include("QuadTree.jl")
include("Radiation.jl")
include("Streamline.jl")
include("Utils.jl")
include("Wind.jl")

end #module