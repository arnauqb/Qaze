export initialize_grids, update_density_and_fm_lines, erase_line_from_grid,
    update_mdot_grid, update_taux_and_xi_grid, refine_density_grid, fill_density_and_fm_grid
using Statistics

function initialize_grids(config::Dict, qsosed)
    disk_r_min = config["disk"]["inner_radius"]
    disk_r_max = config["disk"]["outer_radius"]
    n_r = config["grids"]["n_r"]
    n_z = config["grids"]["n_z"]
    n_disk = config["grids"]["n_r_disk"]
    r_min = config["grids"]["r_min"]
    z_min = config["grids"]["z_min"]
    r_max = config["grids"]["r_max"]
    z_max = config["grids"]["z_max"]
    n_lines = config["wind"]["number_streamlines"]
    if config["grids"]["log_spaced"]
        r_range = logrange(r_min, r_max, n_r)
        z_range = logrange(z_min, z_max, n_z) 
    else
        r_range = range(r_min, stop=r_max, length=n_r)
        z_range = range(z_min, stop=z_max, length=n_z)
    end
    disk_range = 10 .^(range(log10(disk_r_min), stop=log10(disk_r_max), length=n_disk))
    grids = GridsStruct(
        n_r,
        n_z,
        n_disk,
        n_lines,
        r_range,
        z_range,
        sqrt(r_max^2 + z_max^2), # d_max
        disk_range,
        config["wind"]["n_shielding"] * ones(Float64, n_r, n_z), #density
        zeros(Float64, n_lines, n_r, n_z), # density lines
        zeros(Float64, n_r, n_z), #tau_x
        zeros(Float64, n_r, n_z), #ionization
        zeros(Float64, n_r, n_z), #fm
        zeros(Float64, n_lines, n_r, n_z), # fm lines
        config["bh"]["mdot"].* ones(Float64, n_disk), #mdot
        ones(Float64, n_disk), #uv fraction
    )
    return grids
end

function refine_density_grid(r_range_old, r_range_new, z_range_old, z_range_new, density_old)
    n_r = length(r_range_new)
    n_z = length(z_range_new)
    density_new  = zeros(Float64, n_r, n_z)
    for i in 1:n_r
        r = r_range_new[i]
        r_old_arg = searchsortedfirst(r_range_old, r)
        for j in 1:n_z
            z = z_range_new[j]
            z_old_arg = searchsortedfirst(z_range_old, z)
            density_interp = 0.
            counter = 0
            for k1 in 0:2
                if (r_old_arg + k1) > length(r_range_old)
                    continue
                end
                for k2 in -1:1
                    if ((z_old_arg + k2) > length(z_range_old) || (z_old_arg + k2) < 1)
                        continue
                    end
                    density_interp += log10(density_old[r_old_arg + k1, z_old_arg + k2])
                    counter += 1
                end
            end
            density_interp = 10^(density_interp / counter)
            density_new[i,j] = density_interp
        end
    end
    return density_new
end

function refine_density_grid(wind::WindStruct)
    println("Computing optical depths ...")
    delta_taus = zeros(Float64, wind.grids.n_r, wind.grids.n_z)
    for i in 1:wind.grids.n_r - 1
        r1 = wind.grids.r_range[i]
        r2 = wind.grids.r_range[i+1]
        for j in 1:wind.grids.n_z - 1
            z1 = wind.grids.z_range[j]
            z2 = wind.grids.z_range[j+1]
            delta_d = sqrt((r2-r1)^2 + (z2-z1)^2) * wind.bh.R_g
            density = wind.grids.density[i,j]
            delta_tau = delta_d * SIGMA_T * density
            delta_taus[i,j] = delta_tau
        end
    end
    println("refining grid...")
    opt_thick_locations = findall(x -> x>0.1, delta_taus)
    if length(opt_thick_locations) == 0
        println("nothing to refine")
        return nothing
    end
    r_range_new = copy(wind.grids.r_range)
    z_range_new = copy(wind.grids.z_range)
    r_arg_list = []
    z_arg_list = []
    for (i, row) in enumerate(eachrow(opt_thick_locations))
        r_arg, z_arg = Tuple(row[1])
        if i == 1
            push!(r_arg_list, r_arg)
            push!(z_arg_list, z_arg)
        else
            if !(r_arg == r_arg_list[end])
                push!(r_arg_list, r_arg)
            end
            if !(z_arg == z_arg_list[end])
                push!(z_arg_list, z_arg)
            end
        end
    end
    elements_attached = 0
    for r_arg in r_arg_list
        if r_arg == 1
            r = wind.grids.r_range[1]
            rnext = wind.grids.r_range[2]
            r2 = r + (rnext - r) / 2.
            insert!(r_range_new, 2, r2)
            elements_attached += 1
            continue
        end
        r_arg_corr = r_arg + elements_attached
        r = r_range_new[r_arg_corr]
        rprevious = r_range_new[r_arg_corr - 1]
        rnext = r_range_new[r_arg_corr + 1]
        if wind.config["grids"]["log_spaced"]
            r1 = rprevious + 10^((log10(r) - log10(rprevious))/ 2.)
            r2 = r + 10^((log10(rnext) - log10(r))/ 2.)
        else
            r1 = rprevious + (r - rprevious) / 2.
            r2 = r + (rnext - r) / 2.
        end
        deleteat!(r_range_new, r_arg_corr)
        insert!(r_range_new, r_arg_corr, r1)
        insert!(r_range_new, r_arg_corr+1, r2)
        #println("r: $r, $rnext, $r2")
        elements_attached += 1
    end
    elements_attached = 0
    for z_arg in z_arg_list
        if z_arg == 1
            z = wind.grids.z_range[1]
            znext = wind.grids.z_range[2]
            z2 = z + (znext - z) / 2.
            insert!(z_range_new, 2, z2)
            elements_attached += 1
            continue
        end
        z_arg_corr = z_arg + elements_attached
        z = z_range_new[z_arg_corr]
        zprevious = z_range_new[z_arg_corr - 1]
        znext = z_range_new[z_arg_corr + 1]
        if wind.config["grids"]["log_spaced"]
            z1 = zprevious + 10^((log10(z) - log10(zprevious))/ 2.)
            z2 = z + 10^((log10(znext) - log10(z))/ 2.)
        else
            z1 = zprevious + (z - zprevious) / 2.
            z2 = z + (znext - z) / 2.
        end
        deleteat!(z_range_new, z_arg_corr)
        insert!(z_range_new, z_arg_corr, z1)
        insert!(z_range_new, z_arg_corr+1, z2)
        elements_attached += 1
    end
    newsize_r = length(r_range_new)
    newsize_z = length(z_range_new)
    density_new = refine_density_grid(wind.grids.r_range, r_range_new, wind.grids.z_range, z_range_new, wind.grids.density)
    grids = GridsStruct(
        newsize_r,
        newsize_z,
        wind.grids.n_disk,
        wind.grids.n_lines,
        r_range_new,
        z_range_new,
        sqrt(wind.grids.r_range[end]^2 + wind.grids.z_range[end]^2), # d_max
        wind.grids.disk_range,
        density_new, #wind.config["wind"]["n_shielding"] * ones(Float64, newsize_r, newsize_z), #density
        zeros(Float64, wind.grids.n_lines, newsize_r, newsize_z), # density lines
        zeros(Float64, newsize_r, newsize_z), #tau_x
        zeros(Float64, newsize_r, newsize_z), #ionization
        zeros(Float64, newsize_r, newsize_z), #fm
        zeros(Float64, wind.grids.n_lines, newsize_r, newsize_z), # fm lines
        wind.config["bh"]["mdot"].* ones(Float64, wind.grids.n_disk), #mdot
        ones(Float64, wind.grids.n_disk), #uv fraction
    )
    wind.grids = grids
    #println("filling grids...")
    #fill_density_and_fm_grid(wind)
    return nothing
end

function fill_density_and_fm_grid(wind::WindStruct)
    for (line_id, line) in enumerate(wind.lines)
        for i in 1:(size(line.p.u_hist)[1]-1)
            r1, z1 = line.p.u_hist[i, 1:2]
            rho = line.p.n_hist[i]
            fm = line.p.fm_hist[i]
            r2, z2 = line.p.u_hist[i+1, 1:2]
            lw = r1 / line.p.r_0 * line.p.line_width
            update_density_and_fm_lines(r1, r2, z1, z2, lw, rho, fm, line_id, wind)
        end
    end
end

function update_density_and_fm_lines(r1, r2, z1, z2, lw, rho, fm, line_id, wind::WindStruct)
    if z1 > z2
        z1, z2 = z2, z1
    end
    if r1 > r2
        r1, r2 = r2, r1
    end
    z1_arg = get_index(wind.grids.z_range, z1)
    z2_arg = get_index(wind.grids.z_range, z2)
    if z1_arg == z2_arg
        rmin_arg = get_index(wind.grids.r_range, r1)
        rmax_arg = get_index(wind.grids.r_range, r2)
        wind.grids.density_lines[line_id, rmin_arg:rmax_arg, z1_arg] .= rho
        wind.grids.fm_lines[line_id, rmin_arg:rmax_arg, z1_arg] .= fm
    else
        m = (r2 - r1) / (z2 - z1)
        n = r2 - m * z2
        for z_arg in z1_arg:z2_arg
            z = wind.grids.z_range[z_arg]
            r = m * z + n
            rmin_arg = get_index(wind.grids.r_range, r-lw)
            rmax_arg = get_index(wind.grids.r_range, r+lw)
            wind.grids.density_lines[line_id, rmin_arg:rmax_arg, z_arg] .= rho
            wind.grids.fm_lines[line_id, rmin_arg:rmax_arg, z_arg] .= fm
            for r_arg in rmin_arg:rmax_arg
                idx_nz = findall(!iszero, wind.grids.density_lines[:, r_arg, z_arg])
                if length(idx_nz) == 0
                    continue
                else
                    values = wind.grids.density_lines[idx_nz, r_arg, z_arg]
                    rho_mean = mean(values)
                    wind.grids.density[r_arg, z_arg] = rho_mean
                end
                idx_nz = findall(!iszero, wind.grids.fm_lines[:, r_arg, z_arg])
                if length(idx_nz) == 0
                    continue
                else
                    values = wind.grids.fm_lines[idx_nz, r_arg, z_arg]
                    fm_mean = mean(values)
                    wind.grids.fm[r_arg, z_arg] = fm_mean
                end
            end
        end
    end
end

function erase_line_from_grid(line_id, wind::WindStruct)
    wind.grids.density_lines[line_id, :, :] .= 0.
    wind.grids.fm_lines[line_id, :, :] .= 0.
    for i in 1:wind.grids.n_r
        for j in wind.grids.n_z 
            wind.grids.density[i,j] = mean(wind.grids.density_lines[:, i, j]) 
            wind.grids.fm[i,j] = mean(wind.grids.fm_lines[:, i, j])
        end
    end
end

function update_mdot_grid(wind::WindStruct)
    accumulated_wind = 0.
    for line in reverse(wind.lines)
        r_0 = line.p.r_0
        width = line.p.line_width
        rmin_arg = get_index(wind.grids.disk_range, r_0 - width/2.)
        rmax_arg = get_index(wind.grids.disk_range, r_0 + width/2.)
        mw = compute_line_mdot(line, wind)
        mw_norm = mw / mass_accretion_rate(wind.bh)
        accumulated_wind += mw_norm
        mdot = wind.bh.mdot - accumulated_wind
        wind.grids.mdot[rmin_arg:rmax_arg] .= max(mdot , 0.)
    end
    rmin = wind.config["wind"]["initial_radius"]
    rmin_arg = get_index(wind.grids.disk_range, rmin )
    wind.grids.mdot[1:rmin_arg] .= max(wind.bh.mdot - accumulated_wind, 0.)
    return wind.grids.mdot
end

function update_taux_and_xi_grid(wind::WindStruct)
    #rz_grid = collect(Iterators.product(1:wind.grids.n_r, 1:wind.grids.n_z)) 
    #f(x) = compute_tau_x(wind.grids.r_range[x[1]], wind.grids.z_range[x[2]], wind)
    #println(f((1, 10)))
    #wind.grids.tau_x[:,:] = f.(rz_grid)
    #f(x) = ionization_parameter(wind.grids.r_range[x[1]], wind.grids.z_range[x[2]],
    #                            wind.grids.density[x[1],x[2]], wind.grids.tau_x[x[1],x[2]], wind)
    #wind.grids.ionization[:,:] = f.(rz_grid)
    for i in 1:wind.grids.n_r
        r = wind.grids.r_range[i]
        for j in 1:wind.grids.n_z
            z = wind.grids.z_range[j]
            taux = compute_tau_x(r, z, wind)
            density = wind.grids.density[i,j]
            xi = ionization_parameter(r, z, density, taux, wind)
            wind.grids.tau_x[i, j] = taux
            wind.grids.ionization[i, j] = xi
        end
    end
    return nothing
end

