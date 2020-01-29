include("constants.jl")
include("structures.jl")

function gravity(r, z, bh::BlackHole)
    #=
    Gravitational force in units of Rg/c^2
    =#
    d = sqrt(r^2 + z^2)
    constant = -1. / d^2
    result = constant * [r / d, z / d]
    return result
end




