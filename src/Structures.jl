using PyCall
using RegionTrees: AbstractRefinery, Cell
using Interpolations
using VoronoiDelaunay
import VoronoiDelaunay: getx, gety
export GridsStruct, BlackHoleStruct, RadiationStruct, StreamlineStruct, WindStruct, CellData, CustomPoint

struct GridsStruct
    n_disk::Int64
    n_lines::Int64
    r_min::Float64
    r_max::Float64
    z_min::Float64
    z_max::Float64
    disk_range::Array{Float64,1}
    mdot::Array{Float64,1}
    uv_fractions::Array{Float64,1}
    n_vacuum::Float64
end

struct BlackHoleStruct
    M::Float64
    mdot::Float64
    spin::Float64
    isco::Float64
    eta:: Float64
    R_g::Float64
    disk_r_min::Float64
    disk_r_max::Float64
end

mutable struct RadiationStruct
    bol_luminosity::Float64
    eddington_luminosity::Float64
    f_uv::Float64
    f_x::Float64
    xray_luminosity::Float64
    force_constant::Float64
    eta_interpolator::Interpolations.Extrapolation
    k_interpolator::Interpolations.Extrapolation
    include_tauuv::Bool
    include_fm_tauuv::Bool
end

struct CustomPoint <: AbstractPoint2D
    _x::Float64
    _y::Float64
    _density::Float64
    _linewidth::Float64
    _type::String
    _m::Float64
    _n::Float64
end
CustomPoint(x::Float64, y::Float64) = CustomPoint(x, y, 1e2, 0.0, "out", 0.0, 0.0)
getx(p::CustomPoint) = p._x
gety(p::CustomPoint) = p._y

mutable struct WindStruct
    config::Dict
    bh::BlackHoleStruct
    sed::PyObject
    grids::GridsStruct
    quadtree::Cell
    tessellation::DelaunayTessellation2D{CustomPoint}
    radiation::RadiationStruct
    lines::Array{Any,1}
    interpolators::Array{Any, 1}
    lines_initial_radius::Float64
    lines_range::Array{Float64,1}
    lines_widths::Array{Float64,1}
    z_0::Float64
    r_max::Float64
    z_max::Float64
end

mutable struct StreamlineStruct
    wind::WindStruct
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
    outofdomain::Bool
    crossing_counter::Int64
    u_hist::Array{Float64,2}
    n_hist::Array{Float64,1}
    tau_x_hist::Array{Float64,1}
    fm_hist::Array{Float64,1}
    xi_hist::Array{Float64,1}
    dv_dr_hist::Array{Float64,1}
    a_r_hist::Array{Float64,1}
    a_z_hist::Array{Float64,1}
end

struct CellData
    line_id::Array{Int,1}
    z_positions::Array{Float64, 1}
    densities::Array{Float64, 1}
    taus::Array{Float64, 1}
    #z_max::Float64
    #direction::Int # 1 up, -1 down, 0 undecided
end
