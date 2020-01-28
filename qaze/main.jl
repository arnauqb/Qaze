using ConfParser
include("constants.jl")
include("black_hole.jl")
include("grid.jl")
struct Wind
    bh::BlackHole
    grids::Grids
end
function readconfig(config_file::String)
    config = Dict()
    conf = ConfParse(config_file)
    parse_conf!(conf)
    config["general"] = retrieve(conf, "general")
    config["bh"] = retrieve(conf, "bh")
    config["disk"] = retrieve(conf, "disk")
    config["grids"] = retrieve(conf, "grids")
    config["radiation"] = retrieve(conf, "radiation")
    config["wind"] = retrieve(conf, "wind")
    return config
end

function initialize(config::Dict)
    M = parse.(Float64, config["bh"]["M"][1])
    mdot = parse(Float64, config["bh"]["mdot"][1])
    spin = parse(Float64, config["bh"]["spin"][1])
    eta = parse(Float64, config["bh"]["eta"][1])
    R_g = G * M * M_SUN / (C^2)
    bh = BlackHole(M, mdot, spin, eta, R_g)
    n_r = parse(Int64, config["grids"]["n_r"][1])
    n_z = parse(Int64, config["grids"]["n_z"][1])
    n_disk = parse(Int64, config["grids"]["n_r_disk"][1])
    r_min = parse(Float64, config["grids"]["r_min"][1])
    z_min = parse(Float64, config["grids"]["z_min"][1])
    r_max = parse(Float64, config["grids"]["r_max"][1])
    z_max = parse(Float64, config["grids"]["z_max"][1])
    disk_r_min = parse(Float64, config["disk"]["inner_radius"][1])
    disk_r_max = parse(Float64, config["disk"]["outer_radius"][1])
    r_range = 10 .^(range(log10(r_min), stop=log10(r_max), length=n_r))
    z_range = 10 .^(range(log10(z_min), stop=log10(z_max), length=n_z))
    disk_range = 10 .^(range(log10(disk_r_min), stop=log10(disk_r_max), length=n_disk))
    grids = Grids(
        n_r,
        n_z,
        n_disk,
        r_range,
        z_range,
        disk_range,
        zeros(Float64, n_r, n_z), #density
        zeros(Float64, n_r, n_z), #tau_x
        zeros(Float64, n_r, n_z), #ionization
    )
    wind = Wind(bh, grids)
    return wind
end

config = readconfig("config.ini")
wind = initialize(config)

