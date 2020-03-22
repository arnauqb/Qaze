export initialize
using PyCall
using RegionTrees
using Interpolations


"Initialize wind model from config file in TOML format"
function initialize(config::String)
    config = TOML.parsefile(config)
    initialize(config)
end

"Initialize wind model from Dictionary"
function initialize(config::Dict)
    bh = initialize_blackhole(config)
    sed = initialize_qsosed(config)
    grids = initialize_grids(config, sed)
    radiation = initialize_radiation(config, bh, sed)
    quadtree = initialize_quadtree(config, bh, grids)
    wind = initialize_wind(config, bh, sed,grids, quadtree, radiation)
    initialize_uv_fraction!(wind)
    initialize_json(config["general"]["save_path"], wind)
    return wind
end

function initialize_blackhole(config::Dict)
    M = config["bh"]["M"] 
    mdot =config["bh"]["mdot"]
    spin = config["bh"]["spin"]
    eta = config["bh"]["eta"]
    R_g = G * M * M_SUN / (C^2)
    disk_r_min = config["disk"]["inner_radius"]
    disk_r_max = config["disk"]["outer_radius"]
    bh = BlackHoleStruct(M, mdot, spin, 6., eta, R_g, disk_r_min, disk_r_max)
    return bh
end

function initialize_qsosed(config::Dict)
    qsosed = pyimport("qsosed")
    sed = qsosed.SED(
        M=config["bh"]["M"], 
        mdot=config["bh"]["mdot"],
        number_bins_fractions=config["grids"]["n_r_disk"],
        )
    return sed
end

function initialize_grids(config::Dict, qsosed)
    disk_r_min = config["disk"]["inner_radius"]
    disk_r_max = config["disk"]["outer_radius"]
    n_disk = config["grids"]["n_r_disk"]
    r_min = config["grids"]["r_min"]
    z_min = config["grids"]["z_min"]
    r_max = config["grids"]["r_max"]
    z_max = config["grids"]["z_max"]
    n_vacuum = config["grids"]["n_vacuum"]
    n_lines = config["wind"]["number_streamlines"]
    disk_range = 10 .^(range(log10(disk_r_min), stop=log10(disk_r_max), length=n_disk))
    grids = GridsStruct(
        n_disk,
        n_lines,
        r_min,
        r_max,
        z_min,
        z_max,
        disk_range,
        config["bh"]["mdot"].* ones(Float64, n_disk), #mdot
        ones(Float64, n_disk), #uv fraction
        n_vacuum,
        )
    return grids
end

function initialize_radiation(config::Dict, bh::BlackHoleStruct, sed::PyObject)
    edd_lumin = eddington_luminosity(bh)
    f_uv = config["radiation"]["f_uv"]
    f_x = config["radiation"]["f_x"]
    if (f_uv == "auto" || f_x == "auto")
        f_uv = sed.uv_fraction
        f_x = sed.xray_fraction
    end
    bol_lumin = bh.mdot * edd_lumin
    xray_lumin = f_x * bol_lumin
    force_constant = 3 / (8 * pi * bh.eta)
    fm_k_interpolator = extrapolate(interpolate((K_INTERP_XI_VALUES,), K_INTERP_K_VALUES, Gridded(Linear())), Flat())
    fm_eta_interpolator = extrapolate(interpolate((ETAMAX_INTERP_XI_VALUES,), ETAMAX_INTERP_ETAMAX_VALUES, Gridded(Linear())), Flat())
    rad = RadiationStruct(bol_lumin, edd_lumin, f_uv, f_x, xray_lumin, force_constant, fm_eta_interpolator, fm_k_interpolator, false, false)
    return rad
end

function initialize_quadtree(config::Dict, bh::BlackHoleStruct, grids::GridsStruct)
    quadtree_max_radius = config["grids"]["r_max"]
    quadtree_max_height = config["grids"]["z_max"]
    #quadtree = Cell(SVector(0., 0.), 
    #                SVector(2 * quadtree_max_radius, 2* quadtree_max_height),
    #                [grids.n_vacuum * quadtree_max_height, quadtree_max_height] # density, fm , cell optical thickness, lines that pass through this cell
    #                )
    quadtree = Cell(SVector(0., 0.), 
                    SVector(2 * quadtree_max_radius, 2* quadtree_max_height),
                    CellData(Int64[], Array{Float64, 2}(undef, 2, 0), Float64[], 0, 0)
                    #[0, Array{Float64}(undef, 2, 0), Float64[]] # line_id position density 
                    )
    return quadtree
end

function initialize_wind(config::Dict, bh::BlackHoleStruct, sed::PyObject, grids::GridsStruct, 
                         quadtree::Cell, radiation::RadiationStruct)
    n_lines = config["wind"]["number_streamlines"]
    lines_initial_radius = config["wind"]["initial_radius"]
    if lines_initial_radius == "warm_radius"
        lines_initial_radius = sed.warm_radius
    end
    lines_final_radius = config["wind"]["final_radius"]
    if config["wind"]["log_spaced"]
        line_delimiters = 10 .^ range(log10(lines_initial_radius), log10(lines_final_radius), length=n_lines+1)
        lines_range = []
        for i in 1:n_lines
            r0 = line_delimiters[i] + (line_delimiters[i+1] - line_delimiters[i]) / 2.
            push!(lines_range, r0)
        end
        lines_widths = diff(line_delimiters)
        #dlogr = (log10(lines_final_radius) - log10(lines_initial_radius)) / n_lines
        #log_r_range = [log10(lines_initial_radius) + (i+0.5) * dlogr for i in 0:n_lines-1]
        #lines_range = 10 .^ log_r_range
        #border_value = 10^(log10(lines_range[end]) + dlogr)
        #lines_widths = diff([lines_range ; border_value])
    else
        dr = (lines_final_radius - lines_initial_radius) / n_lines 
        lines_range = [lines_initial_radius + (i+0.5) * dr for i in 0:n_lines-1]
        lines_widths = diff([lines_range ; lines_range[end]])
    end
    lines = Array{Any,1}(undef, config["wind"]["number_streamlines"])
    interpolators = Array{Any,1}(undef, config["wind"]["number_streamlines"])
    z_0 = config["wind"]["z_0"]
    wind = WindStruct(config,
                      bh,
                      sed,
                      grids,
                      quadtree,
                      radiation,
                      lines,
                      interpolators,
                      lines_initial_radius,
                      lines_range,
                      lines_widths,
                      z_0,
                      )
    return wind
end
