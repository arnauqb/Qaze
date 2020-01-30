#using PyCall
struct Grids
    n_r::Int64
    n_z::Int64
    n_disk::Int64
    n_lines::Int64
    r_range::Array{Float64,1}
    z_range::Array{Float64,1}
    d_max::Float64
    disk_range::Array{Float64,1}
    density::Array{Float64,2}
    density_lines::Array{Float64, 3}
    tau_x::Array{Float64,2}
    ionization::Array{Float64,2}
    fm::Array{Float64,2}
    fm_lines::Array{Float64, 3}
    mdot::Array{Float64,1}
    uv_fractions::Array{Float64,1}
end
struct BlackHole
    M::Float64
    mdot::Float64
    spin::Float64
    isco::Float64
    eta:: Float64
    R_g::Float64
    disk_r_min::Float64
    disk_r_max::Float64
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
    sed::PyObject
    grids::Grids
    radiation::Radiation
    lines::Array{Any,1}
    lines_range::Array{Float64,1}
    lines_widths::Array{Float64,1}
end

mutable struct Streamline
    wind::Wind
    line_id::Int64
    r_0::Float64
    z_0::Float64
    v_0::Float64
    v_phi_0::Float64
    n_0::Float64
    v_th::Float64
    l::Float64 # specific angular momentum
    line_width::Float64
    escaped::Bool
    crossing_counter::Int64
    is_first_iter::Bool
    u_hist::Array{Float64,2}
    n_hist::Array{Float64,1}
    fm_hist::Array{Float64,1}
    xi_hist::Array{Float64,1}
end