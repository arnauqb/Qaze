include("constants.jl")
include("structures.jl")

function gravity(bh::BlackHole, r, z, height=0.)
    d = sqrt(r^2 + z^2)
    constant = G * bh.M / d^2
    result = constant * [r / d, z / d]
    return result
end




