struct Grids
    n_r::Int64
    n_z::Int64
    n_disk::Int64
    r_range::Array{Float64,1}
    z_range::Array{Float64,1}
    disk_range::Array{Float64,1}
    density::Array{Float64,2}
    tau_x::Array{Float64,2}
    ionization::Array{Float64,2}
end


