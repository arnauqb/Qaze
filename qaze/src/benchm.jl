using BenchmarkTools
include("radiation.jl")
include("qaze.jl")

r_range = range(6, stop=1600, length=50)
z_range = range(1, stop=500, length=50)
integ = zeros(Float64, length(r_range), length(z_range),2)


function run_int()
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            integ[i,j, 1] = integrate_r(r, z, wind, include_tau_uv=true)
            integ[i,j, 2] = integrate_z(r, z, wind, include_tau_uv=true)
        end
    end
end

function run_int2()
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            integ[i,j,:] = integrate(r, z, wind)
        end
    end
end


#@benchmark run_int