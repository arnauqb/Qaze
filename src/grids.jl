export update_density_and_fm_lines
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