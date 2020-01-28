struct Grids
    n_r::Int64
    n_z::Int64
    n_disk::Int64
    r_range::Array{Float64,1}
    z_range::Array{Float64,1}
    d_max::Float64
    disk_range::Array{Float64,1}
    density::Array{Float64,2}
    tau_x::Array{Float64,2}
    ionization::Array{Float64,2}
    fm::Array{Float64,2}
    mdot::Array{Float64,1}
    uv_fraction::Array{Float64,1}
end
struct BlackHole
    M::Float64
    mdot::Float64
    spin::Float64
    isco::Float64
    eta:: Float64
    R_g::Float64
end

struct Radiation
    bol_luminosity::Float64
    eddington_luminosity::Float64
    f_uv::Float64
    f_x::Float64
    xray_luminosity::Float64
    force_constant::Float64
end
struct Wind
    config::Dict
    bh::BlackHole
    grids::Grids
    radiation::Radiation
    is_first_iter::Bool
end

mutable struct Streamline
    wind::Wind
    r_0::Float64
    z_0::Float64
    v_0::Float64
    v_phi_0::Float64
    n_0::Float64
    v_th::Float64
    l::Float64 # specific angular momentum
    escaped::Bool
end

