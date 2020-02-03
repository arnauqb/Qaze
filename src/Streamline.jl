export compute_density, initialize_line, compute_initial_acceleration, residual!, condition, affect, save, solve_line!
using DifferentialEquations
using Sundials

function compute_density(r, z, v_T, line::StreamlineStruct)
    d = sqrt(r^2 + z^2)
    radial = (line.r_0 / d)^2
    v_ratio = line.v_0 / v_T
    n = line.n_0 * radial * v_ratio
    return n
end

function initialize_line(line_id, r_0, z_0, v_0, n_0, v_th, wind::WindStruct, is_first_iter)
    v_phi_0 = sqrt(1. / r_0)
    l = v_phi_0 * r_0
    a_r_0, a_z_0, fm, xi, dv_dr, tau_x = compute_initial_acceleration(r_0, z_0, v_0, n_0, v_th, l, wind, is_first_iter)
    u0 = [r_0, z_0, 0., v_0]
    du0 = [0., v_0, a_r_0, a_z_0]
    lw = wind.lines_widths[line_id]
    u_hist = reshape(u0, (1,4))
    line = StreamlineStruct(wind, line_id, r_0, z_0, v_0, v_phi_0, n_0,
                            v_th, l, lw, false, 0, is_first_iter, u_hist,
                            [n_0], [tau_x], [fm], [xi], [dv_dr], [a_r_0], [a_z_0])
    tspan = (0., 1e8)
    termination_cb = DiscreteCallback(condition, affect, save_positions=(false, false))
    saved_values_type = SavedValues(Float64, Array{Float64,1})
    saving_cb = SavingCallback(save, saved_values_type)
    cb = CallbackSet(termination_cb, saving_cb)
    problem = DAEProblem(residual!, du0, u0, tspan, p=line, differential_vars=[true, true, true, true])
    integrator = init(problem, IDA(), callback=cb)
    integrator.opts.abstol = 0.
    integrator.opts.reltol = wind.config["wind"]["solver_rtol"]
    return integrator
end

function solve_line!(line)
    solve!(line)
end

function compute_initial_acceleration(r_0, z_0, v_0, n_0, v_th, l, wind::WindStruct, is_first_iter)
    fg = gravity(r_0, z_0, wind.bh)
    tau_x = compute_tau_x(r_0, z_0, wind)
    xi = ionization_parameter(r_0, z_0, n_0, tau_x, wind)
    fm = force_multiplier(1, xi)
    fr = force_radiation(r_0, z_0, fm, wind, include_tau_uv=!is_first_iter)
    centrifugal_term = l^2 / r_0^3
    a_r = fg[1] + fr[1] + centrifugal_term
    a_z = fg[2] + fr[2] 
    #second estimation
    a_T = sqrt(a_r^2 + a_z^2)
    dv_dr = a_T / v_0
    tau_eff = compute_tau_eff(n_0, dv_dr, v_th)
    fm = force_multiplier(tau_eff, xi)
    fr = force_radiation(r_0, z_0, fm, wind, include_tau_uv=!is_first_iter)
    a_r = fg[1] + fr[1] + centrifugal_term
    a_z = fg[2] + fr[2] 
    return [a_r, a_z, fm, xi, dv_dr, tau_x]
end

function residual!(resid, du, u, line, t)
    r, z, v_r, v_z = u
    r_dot, z_dot, v_r_dot, v_z_dot = du
    fg = gravity(r, z, line.wind.bh)
    v_T = sqrt(r_dot^2 + z_dot^2)
    a_T = sqrt(v_r_dot^2 + v_z_dot^2)
    dv_dr = a_T / v_T
    n = compute_density(r, z, v_T, line)
    tau_x = compute_tau_x(r, z, line.wind)
    xi = ionization_parameter(r, z, n, tau_x, line.wind)
    tau_eff = compute_tau_eff(n, dv_dr, line.v_th)
    fm = force_multiplier(tau_eff, xi)
    fr = force_radiation(r, z, fm, line.wind, include_tau_uv=!line.is_first_iter)
    centrifugal_term = line.l^2 / r^3
    a_r = fg[1] + fr[1] + centrifugal_term
    a_z = fg[2] + fr[2] 
    resid[1] = r_dot - v_r
    resid[2] = z_dot - v_z
    resid[3] = v_r_dot - a_r
    resid[4] = v_z_dot - a_z
end

function condition(u, t, integrator)
    r, z, v_r, v_z = u
    _, _, a_r, a_z = integrator.du
    d = sqrt(r^2 + z^2)
    v_T = sqrt(v_r^2 + v_z^2)
    v_esc = sqrt(2. / d)
    if v_T > v_esc
        integrator.p.escaped = true
    end
    crossing_condition = false
    if r < integrator.p.r_0
        integrator.p.crossing_counter += 1
        if integrator.p.crossing_counter > 4
            crossing_condition = true
        end
    end
    
    escaped_condition = d > integrator.p.wind.grids.d_max
    failed_condtion = z < integrator.p.z_0 
    cond = escaped_condition | failed_condtion | crossing_condition
    return cond
end

function affect(integrator)
    if integrator.p.escaped
        print(" \U1F4A8")
    else
        print(" \U1F4A5")
    end
    terminate!(integrator)
end

function save(u, t, integrator)
    r, z, v_r, v_z = u
    du = get_du(integrator)
    v_r_dot, v_z_dot, a_r, a_z = du
    v_T = sqrt(v_r^2 + v_z^2)
    a_T = sqrt(a_r^2 + a_z^2)
    dv_dr = a_T / v_T
    n = compute_density(r, z, v_T, integrator.p)
    tau_x = compute_tau_x(r, z, integrator.p.wind)
    xi = ionization_parameter(r, z, n, tau_x, integrator.p.wind)
    tau_eff = compute_tau_eff(n, dv_dr, integrator.p.v_th)
    fm = force_multiplier(tau_eff, xi)
    r_0, z_0, v_r_0, v_z_0 = integrator.p.u_hist[end,:]
    if integrator.p.wind.config["radiation"]["tau_uv_include_fm"]
        update_density_and_fm_lines(r_0, r, z_0, z, integrator.p.line_width, n, fm, integrator.p.line_id, integrator.p.wind)
    else
        update_density_and_fm_lines(r_0, r, z_0, z, integrator.p.line_width, n, 0., integrator.p.line_id, integrator.p.wind)
    end
    integrator.p.u_hist = [integrator.p.u_hist ; transpose(u)]
    push!(integrator.p.fm_hist, fm)
    push!(integrator.p.n_hist, n)
    push!(integrator.p.tau_x_hist, tau_x)
    push!(integrator.p.xi_hist, xi)
    push!(integrator.p.dv_dr_hist, dv_dr)
    push!(integrator.p.a_r_hist, a_r)
    push!(integrator.p.a_z_hist, a_z)
    return u
end
#end