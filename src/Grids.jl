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
    opt_thick = findall(x -> x>0.1, delta_taus)
    if length(opt_thick) == 0
        println("nothing to refine")
        return nothing
    end
    r_min_arg, z_min_arg = Tuple(opt_thick[1])
    r_max_arg, z_max_arg = Tuple(opt_thick[end])
    z_max = wind.grids.z_range[z_max_arg]
    z_min = wind.grids.z_range[z_min_arg]
    r_max = wind.grids.r_range[r_max_arg]
    r_min = wind.grids.r_range[r_min_arg]
    # r part
    nr = r_max_arg - r_min_arg 
    r_range_old = wind.grids.r_range
    r_refined_0 = logrange(r_min, r_max, 2 * (nr + 1) )
    r_range = [wind.grids.r_range[1:r_min_arg-1] ; r_refined_0 ; wind.grids.r_range[r_max_arg+1:end]]
    newsize_r = length(r_range)
    # z part
    nz = z_max_arg - z_min_arg 
    z_range_old = wind.grids.z_range
    z_refined_0 = logrange(z_min, z_max, 2 * (nz + 1) )
    z_range = [wind.grids.z_range[1:z_min_arg-1] ; z_refined_0 ; wind.grids.z_range[z_max_arg+1:end]]
    newsize_z = length(z_range)
    grids = GridsStruct(
        newsize_r,
        newsize_z,
        wind.grids.n_disk,
        wind.grids.n_lines,
        r_range,
        z_range,
        sqrt(wind.grids.r_range[end]^2 + wind.grids.z_range[end]^2), # d_max
        wind.grids.disk_range,
        wind.config["wind"]["n_shielding"] * ones(Float64, newsize_r, newsize_z), #density
        zeros(Float64, wind.grids.n_lines, newsize_r, newsize_z), # density lines
        zeros(Float64, newsize_r, newsize_z), #tau_x
        zeros(Float64, newsize_r, newsize_z), #ionization
        zeros(Float64, newsize_r, newsize_z), #fm
        zeros(Float64, wind.grids.n_lines, newsize_r, newsize_z), # fm lines
        wind.config["bh"]["mdot"].* ones(Float64, wind.grids.n_disk), #mdot
        ones(Float64, wind.grids.n_disk), #uv fraction
    )
    wind.grids = grids
    println("filling grids...")
    fill_density_and_fm_grid(wind)
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

