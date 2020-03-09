using Printf

export initialize_line!, start_lines!, compute_line_mdot, compute_wind_mdot,
       compute_kinetic_luminosity, compute_maximum_velocity, start_iteration!, start_line!

function initialize_line!(i, r, wind::WindStruct)
    T = wind.sed.disk_nt_temperature4(r)^(0.25)
    v_th = thermal_velocity(T)
    if wind.config["wind"]["v_0"] == "thermal"
        v_0 = v_th
    else
        v_0 = wind.config["wind"]["v_0"] / C
    end
    n_0 = wind.config["wind"]["n_0"]
    if n_0 == "cak"
        n_0 = cak_density(r, wind)
    end
    z_0 = wind.z_0
    line_integ = initialize_line!(i, r, z_0, v_0, n_0, v_th, wind)
    return line_integ
end

function start_line!(line_id, r_0, wind::WindStruct)
    line = initialize_line!(line_id, r_0, wind)
    wind.lines[line_id] = line
    @printf("\nSolving line %02d of %02d ", line_id, length(wind.lines))
    flush(stdout)
    solve!(line)
    flush(stdout)
    return line
end

function start_iteration!(it_num, until_line, wind::WindStruct)
    if until_line === nothing
        until_line = wind.config["wind"]["number_streamlines"]
    end
    for (i, r_0) in enumerate(wind.lines_range)
        if i > until_line
            break
        end
        if it_num > 1
            erase_line_from_tree!(i, wind)
        end
        line = start_line!(i, r_0, wind)
        write_line(wind.config["general"]["save_path"], line.p, it_num)
    end
end

function start_lines!(wind::WindStruct, until_line = nothing)
    wind.radiation.include_tauuv = false
    for it_num in 1:wind.config["wind"]["iterations"]
        @printf("Iteration %02d of %02d\n", it_num, wind.config["wind"]["iterations"])
        flush(stdout)
        if it_num > 1
            update_mdot_grid!(wind)
            wind.radiation.include_tauuv = true
            wind.config["radiation"]["tau_uv_include_fm"] && (wind.radiation.include_fm_tauuv = true)
        end
        start_iteration!(it_num, until_line, wind)
        #refine_all(wind)
        write_properties_and_grids(wind.config["general"]["save_path"], wind, it_num)
    end
end

function compute_line_mdot(line, wind::WindStruct)
    if line.p.escaped
        n_0 = line.p.n_0
        v_0 = line.p.v_0
        mw = 2pi * line.p.r_0 * line.p.line_width * wind.bh.R_g^2 * n_0 * v_0 * C * M_P
        return mw
    else
        return 0.
    end
end

function compute_wind_mdot(wind::WindStruct)
    mdot_wind = 0.
    for line in wind.lines
        mdot_wind += compute_line_mdot(line, wind)
    end
    return mdot_wind
end

function compute_kinetic_luminosity(wind::WindStruct)
    kin_lumin = 0.
    for line in wind.lines
        v_r_f = line.p.u_hist[end,3]
        v_z_f = line.p.u_hist[end,4]
        v_f = sqrt(v_r_f^2 + v_z_f^2) * C
        kin_lumin += 0.5 * compute_line_mdot(line, wind) * v_f^2
    end
    return kin_lumin
end

function compute_maximum_velocity(wind::WindStruct)
    maxv = 0.
    for line in wind.lines
        if !line.p.escaped
            continue
        end
        v_r_f = line.p.u_hist[end,3]
        v_z_f = line.p.u_hist[end,4]
        v_f = sqrt(v_r_f^2 + v_z_f^2)
        if v_f > maxv
            maxv = v_f
        end
    end
    return maxv
end