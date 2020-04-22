using Printf
using JLD2

export initialize_line!, start_lines!, compute_line_mdot, compute_wind_mdot,
       compute_kinetic_luminosity, compute_maximum_velocity, start_iteration!, start_line!

function initialize_line!(i::Int, wind::WindStruct)
    r = wind.lines_range[i]
    return initialize_line!(i , r, wind)
end

function initialize_line!(i, r, wind::WindStruct ; nolaunch=false)
    if nolaunch
        line_integ = initialize_line!(i, r, 0, 1e-5, 1e8, 1e-5, wind)
        return line_integ
    end
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

function start_iteration!(it_num, wind::WindStruct, from_line = 1, until_line=nothing)
    if until_line === nothing
        until_line = wind.config["wind"]["number_streamlines"]
    end
    if it_num == 1
        wind.radiation.include_tauuv = false 
    elseif it_num == 2
        wind.radiation.include_tauuv = true
    else
        wind.radiation.include_tauuv = true
        update_mdot_grid!(wind)
    end
    for (i, r_0) in enumerate(wind.lines_range)
        if i < from_line
            continue
        end
        if i > until_line
            break
        end
        r0_idx = get_index(wind.grids.disk_range, r_0)
        if wind.grids.mdot[r0_idx] <= 0
            wind.lines[i] = initialize_line!(i, r_0, wind, nolaunch=true)#nothing
            continue
        end

        if it_num > 1
            println("erasing line...")
            flush(stdout)
            quadtree_initialize(wind)
            for (j, line) in enumerate(wind.lines)
                if j == i 
                    continue
                end
                println("filling line $j of $(length(wind.lines))")
                quadtree_fill_line(line, j, wind)
            end

            #erase_line_from_tree!(i, wind)
            #quadtree_erase_line(i, wind)
        end
        # check if there is mass to launch the line
        line = start_line!(i, r_0, wind)
        #println("filling and refining line...")
        #quadtree_fill_line(line, i, wind)
        #fill_and_refine_line!(line, i, wind)
        write_line(wind.config["general"]["save_path"], line.p, it_num)
        #if (i+1) % 3 == 0
        #    @save wind.config["save_jld2"] wind.quadtree wind.grids.mdot
        #end
    end
    mdotcopy = copy(wind.grids.mdot)
    update_mdot_grid!(wind)
    write_properties_and_grids(wind.config["general"]["save_path"], wind, it_num)
    wind.grids.mdot .= mdotcopy
    qt = wind.quadtree
    @save wind.config["general"]["save_jld2"] qt
end

function start_lines!(wind::WindStruct ; from_iter=1, from_line=1, until_line=nothing)
    #wind.radiation.include_tauuv = true#false 
    for it_num in from_iter:wind.config["wind"]["iterations"]
        @printf("Iteration %02d of %02d\n", it_num, wind.config["wind"]["iterations"])
        flush(stdout)
        start_iteration!(it_num, wind, from_line, until_line)
        #refine_all(wind)
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
    for i in 1:length(wind.lines)
        isassigned(wind.lines, i) || continue
        line = wind.lines[i]
        mdot_wind += compute_line_mdot(line, wind)
    end
    return mdot_wind
end

function compute_kinetic_luminosity(wind::WindStruct)
    kin_lumin = 0.
    for i in 1:length(wind.lines)
        isassigned(wind.lines, i) || continue
        line = wind.lines[i]
        v_r_f = line.p.u_hist[end,3]
        v_z_f = line.p.u_hist[end,4]
        v_f = sqrt(v_r_f^2 + v_z_f^2) * C
        kin_lumin += 0.5 * compute_line_mdot(line, wind) * v_f^2
    end
    return kin_lumin
end

function compute_maximum_velocity(wind::WindStruct)
    maxv = 0.
    for i in 1:length(wind.lines)
        isassigned(wind.lines, i) || continue
        line = wind.lines[i]
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
