export compute_density, initialize_line, compute_initial_acceleration, residual!, condition, affect, save, solve_line!
using DifferentialEquations
using Sundials
using Statistics
using RegionTrees


"Updates the density of the streamline giving its current position and velocity,
using mass conservation."
function compute_density(r, z, v_T, line::StreamlineStruct)
    #@assert r >= 0
    #@assert z >= 0
    #d = sqrt(r^2 + z^2)
    #radial = (line.r_0 / d)^2
    #v_ratio = line.v_0 / v_T
    #n = line.n_0 * radial * v_ratio
    #return n
    return line.n_0
end

"Initializes the IDA solver given the initial conditions"
function initialize_line!(line_id, r_0, z_0, v_0, n_0, v_th, wind::WindStruct)
    v_phi_0 = sqrt(1. / r_0) # initial kepler velocity
    l = v_phi_0 * r_0 # initial angular momentum
    lw0 = wind.lines_widths[line_id]
    linewidth_normalized = lw0 / r_0
    fill_and_refine!([r_0, z_0], [r_0, z_0], linewidth_normalized, n_0, 0., wind) # fill first point
    a_r_0, a_z_0, fm, xi, dv_dr, tau_x = compute_initial_acceleration(r_0, z_0, v_0, n_0, v_th, l, wind)
    u0 = [r_0, z_0, 0., v_0]
    du0 = [0., v_0, a_r_0, a_z_0]
    u_hist = reshape(u0, (1,4))
    line = StreamlineStruct(wind,
                            line_id,
                            r_0,
                            z_0,
                            v_0,
                            v_phi_0,
                            n_0,
                            v_th,
                            l,
                            lw0,
                            false,
                            false,
                            0,
                            u_hist,
                            [n_0],
                            [tau_x],
                            [fm],
                            [xi],
                            [dv_dr],
                            [a_r_0],
                            [a_z_0])
    tspan = (0., 1e8)
    termination_cb = DiscreteCallback(condition, affect, save_positions=(false, false))
    saved_values_type = SavedValues(Float64, Array{Float64,1})
    saving_cb = SavingCallback(save, saved_values_type)
    steadystate_cb = TerminateSteadyState(1e-8, 1e-6) 
    cb = CallbackSet(termination_cb, saving_cb, steadystate_cb)
    problem = DAEProblem(residual!, du0, u0, tspan, p=line, differential_vars=[true, true, true, true])
    integrator = init(problem, IDA(), callback=cb)#, dtmax=5e-2 * wind.bh.R_g / C)
    integrator.opts.abstol = 0.
    integrator.opts.reltol = wind.config["wind"]["solver_rtol"]
    return integrator
end

function solve_line!(line)
    solve!(line)
end

"Computes the initial acceleration, consistently with the current values of the ionization parameter."
function compute_initial_acceleration(r_0, z_0, v_0, n_0, v_th, l, wind::WindStruct)
    fg = gravity(r_0, z_0, wind.bh)
    tau_x = compute_tau_x(r_0, z_0, wind)
    xi = ionization_parameter(r_0, z_0, n_0, tau_x, wind)
    fm = force_multiplier(1, xi)
    fr = force_radiation(r_0, z_0, fm, wind, include_tau_uv=wind.radiation.include_tauuv)
    centrifugal_term = l^2 / r_0^3
    a_r = fg[1] + fr[1] + centrifugal_term
    a_z = fg[2] + fr[2] 
    #second estimation
    a_T = sqrt(a_r^2 + a_z^2)
    dv_dr = a_T / v_0
    tau_eff = compute_tau_eff(n_0, dv_dr, v_th)
    fm = force_multiplier(tau_eff, xi)
    fr = force_radiation(r_0, z_0, fm, wind, include_tau_uv=wind.radiation.include_tauuv)
    a_r = fg[1] + fr[1] + centrifugal_term
    a_z = fg[2] + fr[2] 
    return [a_r, a_z, fm, xi, dv_dr, tau_x]
end

"Residual function of the implict diffeq system"
function residual!(resid, du, u, line, t)
    r, z, v_r, v_z = u
    r_dot, z_dot, v_r_dot, v_z_dot = du
    if z < 0 || r < 0 # if integrator tries outside domain, switch off radiation
        fg = gravity(r, z, line.wind.bh)
        centrifugal_term = line.l^2 / r^3
        a_r = fg[1] + centrifugal_term
        a_z = fg[2]
        resid[1] = r_dot - v_r
        resid[2] = z_dot - v_z
        resid[3] = v_r_dot - a_r
        resid[4] = v_z_dot - a_z
        return nothing
    end
    fg = gravity(r, z, line.wind.bh)
    v_T = sqrt(r_dot^2 + z_dot^2)
    a_T = sqrt(v_r_dot^2 + v_z_dot^2)
    dv_dr = a_T / v_T
    n = compute_density(r, z, v_T, line)
    tau_x = compute_tau_x(r, z, line.wind)
    xi = ionization_parameter(r, z, n, tau_x, line.wind)
    tau_eff = compute_tau_eff(n, dv_dr, line.v_th)
    fm = force_multiplier(tau_eff, xi)
    fr = force_radiation(r, z, fm, line.wind, include_tau_uv=line.wind.radiation.include_tauuv)
    centrifugal_term = line.l^2 / r^3
    a_r = fg[1] + fr[1] + centrifugal_term
    a_z = fg[2] + fr[2] 
    resid[1] = r_dot - v_r
    resid[2] = z_dot - v_z
    resid[3] = v_r_dot - a_r
    resid[4] = v_z_dot - a_z
end

"Termination condition"
function condition(u, t, integrator)
    r, z, v_r, v_z = u
    _, _, a_r, a_z = integrator.du
    d = sqrt(r^2 + z^2)
    v_T = sqrt(v_r^2 + v_z^2)
    v_esc = sqrt(2. / d)
    v_T > v_esc && (integrator.p.escaped=true)
    crossing_condition = false
    if r < integrator.p.r_0 - 1
        integrator.p.crossing_counter += 1
        if integrator.p.crossing_counter >= 10
            crossing_condition = true
        end
    end

    #stalling_condition = false
    #z_hist = integrator.p.u_hist[:,2]
    #if length(z_hist) > 200
    #    if std(z_hist[end-10:end]) < 0.02
    #        println("stalled")
    #        stalling_condition = true
    #    end
    #end
    escaped_condition = (r >= integrator.p.wind.grids.r_max) || (z >= integrator.p.wind.grids.z_max)
    escaped_condition && (integrator.p.outofdomain=true)
    failed_condtion = z < integrator.p.z_0  || r < 0.
    cond = escaped_condition | failed_condtion | crossing_condition 
    return cond
end

function affect(integrator)
    if integrator.p.escaped
        print(" \U1F4A8")
    elseif integrator.p.outofdomain
        print(" \U2753")
    else
        print(" \U1F4A5")
    end
    terminate!(integrator)
end

"Saves current iteration data"
function save(u, t, integrator)
    r, z, v_r, v_z = u
    if any([r,z] .< 0)
        terminate!(integrator)
        return u
    end
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
    r_0, z_0, v_r_0, v_z_0 = integrator.p.u_hist[end, :]
    linewidth_normalized = integrator.p.line_width / integrator.p.r_0 
    currentpoint = [r, z]
    previouspoint = [r_0, z_0]
    fill_and_refine!(previouspoint, currentpoint, linewidth_normalized, n, fm, integrator.p.wind)
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