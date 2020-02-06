export initialize
#ENV["PYCALL_JL_RUNTIME_PYTHON"]="/home/arnau/Documents/qwind/env/bin/python"
using PyCall


function initialize(config::String)
    config = TOML.parsefile(config)
    initialize(config)
end

function initialize(config::Dict)
    M = config["bh"]["M"] 
    mdot =config["bh"]["mdot"]
    spin = config["bh"]["spin"]
    eta = config["bh"]["eta"]
    R_g = G * M * M_SUN / (C^2)
    disk_r_min = config["disk"]["inner_radius"]
    disk_r_max = config["disk"]["outer_radius"]
    bh = BlackHoleStruct(M, mdot, spin, 6., eta, R_g, disk_r_min, disk_r_max)
    qsosed = pyimport("qsosed")
    sed = qsosed.SED(M=M, mdot=mdot, number_bins_fractions=config["grids"]["n_r_disk"])
    grids = initialize_grids(config, sed)
    n_lines = config["wind"]["number_streamlines"]
    lines_initial_radius = config["wind"]["initial_radius"]
    if lines_initial_radius == "warm_radius"
        lines_initial_radius = sed.warm_radius
    end
    lines_final_radius = config["wind"]["final_radius"]
    if config["wind"]["log_spaced"]
        dlogr = (log10(lines_final_radius) - log10(lines_initial_radius)) / n_lines
        log_r_range = [log10(lines_initial_radius) + (i+0.5) * dlogr for i in 0:n_lines-1]
        lines_range = 10 .^ log_r_range
        border_value = 10^(log10(lines_range[end]) + dlogr)
        lines_widths = diff([lines_range ; border_value])
    else
        dr = (lines_final_radius - lines_initial_radius) / n_lines 
        lines_range = [lines_initial_radius + (i+0.5) * dr for i in 0:n_lines-1]
        lines_widths = diff([lines_range ; lines_range[end]])
    end
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
    rad = RadiationStruct(bol_lumin, edd_lumin, f_uv, f_x, xray_lumin, force_constant)
    lines = Array{Any,1}(undef, config["wind"]["number_streamlines"])
    wind = WindStruct(config, bh, sed, grids, rad, lines, lines_initial_radius, lines_range, lines_widths)
    initialize_uv_fraction(wind)
    initialize_json(config["general"]["save_path"], wind)
    return wind
end