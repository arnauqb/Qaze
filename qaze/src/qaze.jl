module Qaze
ENV["PYCALL_JL_RUNTIME_PYTHON"]="/home/arnau/Documents/qwind/env/bin/python"
using PyCall
qsosed = pyimport("qsosed")
using TOML
println(pyimport("sys").executable)
include("constants.jl")
include("structures.jl")
include("black_hole.jl")
include("radiation.jl")
include("utils.jl")
include("grids.jl")
include("streamline.jl")

function initialize(config::Dict)
    M = config["bh"]["M"] 
    mdot =config["bh"]["mdot"]
    spin = config["bh"]["spin"]
    eta = config["bh"]["eta"]
    R_g = G * M * M_SUN / (C^2)
    disk_r_min = config["disk"]["inner_radius"]
    disk_r_max = config["disk"]["outer_radius"]
    bh = BlackHole(M, mdot, spin, 6., eta, R_g, disk_r_min, disk_r_max)
    n_r = config["grids"]["n_r"]
    n_z = config["grids"]["n_z"]
    n_disk = config["grids"]["n_r_disk"]
    r_min = config["grids"]["r_min"]
    z_min = config["grids"]["z_min"]
    r_max = config["grids"]["r_max"]
    z_max = config["grids"]["z_max"]
    n_lines = config["wind"]["number_streamlines"]
    r_range = 10 .^(range(log10(r_min), stop=log10(r_max), length=n_r))
    z_range = 10 .^(range(log10(z_min), stop=log10(z_max), length=n_z))
    disk_range = 10 .^(range(log10(disk_r_min), stop=log10(disk_r_max), length=n_disk))
    lines_initial_radius = config["wind"]["initial_radius"]
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
    grids = Grids(
        n_r,
        n_z,
        n_disk,
        n_lines,
        r_range,
        z_range,
        sqrt(r_max^2 + z_max^2), # d_max
        disk_range,
        config["wind"]["n_shielding"] * ones(Float64, n_r, n_z), #density
        zeros(Float64, n_lines, n_r, n_z), # fm lines
        zeros(Float64, n_r, n_z), #tau_x
        zeros(Float64, n_r, n_z), #ionization
        zeros(Float64, n_r, n_z), #fm
        zeros(Float64, n_lines, n_r, n_z), # fm lines
        mdot .* ones(Float64, n_disk), #mdot
        ones(Float64, n_disk), #uv fraction
    )
    sed = qsosed.SED(M=M, mdot=mdot, number_bins_fractions=n_disk)
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
    rad = Radiation(bol_lumin, edd_lumin, f_uv, f_x, xray_lumin, force_constant)
    lines = Array{Any,1}(undef, config["wind"]["number_streamlines"])
    wind = Wind(config, bh, sed, grids, rad, lines, lines_range, lines_widths)
    initialize_uv_fraction(wind)
    initialize_json(config["general"]["save_path"], wind)
    return wind
end

function initialize_line(i, r, wind::Wind, is_first_iter)
    T = wind.sed.disk_nt_temperature4(r)^(0.25)
    v_th = thermal_velocity(T)
    if wind.config["wind"]["v_0"] == "thermal"
        v_0 = v_th
    end
    n_0 = wind.config["wind"]["n_0"]
    z_0 = wind.config["wind"]["z_0"]
    line_integ = initialize_line(i, r, z_0, v_0, n_0, v_th, wind, is_first_iter)
    return line_integ
end

#function initialize_lines(wind::Wind)
#    for (i, r) in enumerate(wind.lines_range)
#        wind.lines[i] = initialize_line(i, r, wind)
#end

function start_lines(wind::Wind)
    is_first_iter = false
    for it_num in 1:wind.config["wind"]["iterations"]
        for (i, r) in enumerate(wind.lines_range)
            line = initialize_line(i, r, wind, is_first_iter)
            println("Solving line $i of $(length(wind.lines))")
            solve!(line)
            write_line(wind.config["general"]["save_path"], line.p, it_num)
        end
        write_properties_and_grids(wind.config["general"]["save_path"], wind, it_num)
        is_first_iter = true
    end
end
config = TOML.parsefile("config.toml")
wind = initialize(config)
start_lines(wind)
#T = wind.sed.disk_nt_temperature4(100)^(0.25)
#v_th = thermal_velocity(T)
#println("T: $T, v_th: $v_th")
#line = initialize_line(1, 100, 1, v_th, 2e8, v_th, wind)
#write_line(json_file)

end #module
