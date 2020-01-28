using TOML
include("constants.jl")
include("structures.jl")
include("black_hole.jl")
include("radiation.jl")
include("utils.jl")
include("grids.jl")
function initialize(config::Dict)
    M = config["bh"]["M"]
    mdot =config["bh"]["mdot"]
    spin = config["bh"]["spin"]
    eta = config["bh"]["eta"]
    R_g = G * M * M_SUN / (C^2)
    bh = BlackHole(M, mdot, spin, 6., eta, R_g)
    n_r = config["grids"]["n_r"]
    n_z = config["grids"]["n_z"]
    n_disk = config["grids"]["n_r_disk"]
    r_min = config["grids"]["r_min"]
    z_min = config["grids"]["z_min"]
    r_max = config["grids"]["r_max"]
    z_max = config["grids"]["z_max"]
    disk_r_min = config["disk"]["inner_radius"]
    disk_r_max = config["disk"]["outer_radius"]

    r_range = 10 .^(range(log10(r_min), stop=log10(r_max), length=n_r))
    z_range = 10 .^(range(log10(z_min), stop=log10(z_max), length=n_z))
    disk_range = 10 .^(range(log10(disk_r_min), stop=log10(disk_r_max), length=n_disk))
    grids = Grids(
        n_r,
        n_z,
        n_disk,
        r_range,
        z_range,
        sqrt(r_max^2 + z_max^2),
        disk_range,
        config["wind"]["rho_shielding"] * ones(Float64, n_r, n_z), #density
        zeros(Float64, n_r, n_z), #tau_x
        zeros(Float64, n_r, n_z), #ionization
        zeros(Float64, n_r, n_z), #fm
        ones(Float64, n_disk), #mdot
        ones(Float64, n_disk), #uv fraction
    )
    edd_lumin = eddington_luminosity(bh)
    f_uv = config["radiation"]["f_uv"]
    f_x = config["radiation"]["f_x"]
    bol_lumin = bh.mdot * edd_lumin
    xray_lumin = f_x * bol_lumin
    force_constant = 3 / (8 * pi * bh.eta)
    rad = Radiation(bol_lumin, edd_lumin, f_uv, f_x, xray_lumin, force_constant)
    wind = Wind(config, bh, grids, rad, false)
    return wind
end

config = TOML.parsefile("config.toml")
wind = initialize(config)
lines = Array{Any, wind.config["wind"]["number_streamlines"]}
