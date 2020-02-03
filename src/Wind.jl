export initialize_line, start_lines, compute_line_mdot, compute_wind_mdot,
       compute_kinetic_luminosity, compute_maximum_velocity

function initialize_line(i, r, wind::WindStruct, is_first_iter)
    T = wind.sed.disk_nt_temperature4(r)^(0.25)
    v_th = thermal_velocity(T)
    if wind.config["wind"]["v_0"] == "thermal"
        v_0 = v_th
    else
        v_0 = v_0 / C
    end
    n_0 = wind.config["wind"]["n_0"]
    if n_0 == "cak"
        n_0 = cak_density(r, wind)
    end
    z_0 = wind.config["wind"]["z_0"]
    line_integ = initialize_line(i, r, z_0, v_0, n_0, v_th, wind, is_first_iter)
    return line_integ
end

function start_lines(wind::WindStruct)
    is_first_iter = true
    for it_num in 1:wind.config["wind"]["iterations"]
        if !is_first_iter
            wind.grids.density_lines[:,:,:] .= 0.
            wind.grids.density[:,:] .= wind.config["wind"]["n_shielding"]
            wind.grids.fm_lines[:,:,:] .= 0.
            wind.grids.fm[:,:] .= 0.
        end
        for (i, r) in enumerate(wind.lines_range)
            if !is_first_iter
                erase_line_from_grid(i, wind)
            end
            line = initialize_line(i, r, wind, is_first_iter)
            wind.lines[i] = line
            print("\nSolving line $i of $(length(wind.lines))")
            solve!(line)
            write_line(wind.config["general"]["save_path"], line.p, it_num)
        end
        if wind.config["general"]["consistent_mdot"]
            update_mdot_grid(wind)
        end
        write_properties_and_grids(wind.config["general"]["save_path"], wind, it_num)
        is_first_iter = false
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