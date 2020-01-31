export update_density_and_fm_lines, erase_line_from_grid, compute_mdot_grid, update_taux_and_xi_grid
using Statistics

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