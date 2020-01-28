include("constants.jl")
struct BlackHole
    M::Float64
    mdot::Float64
    spin::Float64
    eta:: Float64
    R_g::Float64
end

function gravity(bh::BlackHole, r, z, height=0.)
    d = sqrt(r^2 + z^2)
    constant = G * bh.M / d^2
    result = constant * [r / d, z / d]
    return result
end




