using DifferentialEquations
using Sundials
include("radiation.jl")
include("structures.jl")
include("black_hole.jl")


function compute_density(r, z, v_T, line::Streamline)
    d = sqrt(r^2 + z^2)
    radial = (line.r_0 / d)^2
    v_ratio = line.v_0 / v_T
    n = line.n_0 * radial * v_ratio
    return n
end

function initialize_line(r_0, z_0, v_0, n_0, v_th, wind::Wind)
    v_phi_0 = sqrt(1. / r_0)
    l = v_phi_0 * r_0
    a_r_0, a_z_0 = compute_initial_acceleration(r_0, z_0, v_0, n_0, v_th, l, wind)
    line = Streamline(wind, r_0, z_0, v_0, v_phi_0, n_0, v_th, l, false)
    u0 = [r_0, z_0, 0., v_0]
    du0 = [0., v_0, a_r_0, a_z_0]
    tspan = (0., 1e8)
    c_callback = DiscreteCallback(condition, affect, save_positions=(false, false))
    problem = DAEProblem(residual!, du0, u0, tspan, p=line)
    integrator = init(problem, IDA(), callback=c_callback)
    integrator.opts.abstol = 0.
    integrator.opts.reltol = wind.config["wind"]["solver_rtol"]
    return integrator
end

function compute_initial_acceleration(r_0, z_0, v_0, n_0, v_th, l, wind::Wind)
    fg = gravity(r_0, z_0, wind.bh)
    xi = ionization_parameter(r_0, z_0, n_0, wind)
    fm = force_multiplier(1, xi)
    fr = force_radiation(r_0, z_0, fm, wind, include_tau_uv=!wind.is_first_iter)
    centrifugal_term = l^2 / r_0^3
    a_r = fg[1] + fr[1] + centrifugal_term
    a_z = fg[2] + fr[2] 
    #second estimation
    a_T = sqrt(a_r^2 + a_z^2)
    dv_dr = a_T / v_0
    tau_eff = compute_tau_eff(n_0, dv_dr, v_th)
    fm = force_multiplier(tau_eff, xi)
    fr = force_radiation(r_0, z_0, fm, wind, include_tau_uv=!wind.is_first_iter)
    a_r = fg[1] + fr[1] + centrifugal_term
    a_z = fg[2] + fr[2] 
    return [a_r, a_z]
end

function residual!(resid, du, u, line, t)
    r, z, v_r, v_z = u
    r_dot, z_dot, v_r_dot, v_z_dot = du
    fg = gravity(r, z, line.wind.bh)
    v_T = sqrt(r_dot^2 + z_dot^2)
    a_T = sqrt(v_r_dot^2 + v_z_dot^2)
    dv_dr = a_T / v_T
    n = compute_density(r, z, v_T, line)
    xi = ionization_parameter(r, z, n, line.wind)
    tau_eff = compute_tau_eff(n, dv_dr, line.v_th)
    fm = force_multiplier(tau_eff, xi)
    fr = force_radiation(r, z, fm, line.wind, include_tau_uv=!line.wind.is_first_iter)
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
    escaped_condition = d > integrator.p.wind.grids.d_max
    failed_condtion = z < integrator.p.z_0 
    println("$d, $(integrator.p.wind.grids.d_max)")
    cond = escaped_condition | failed_condtion 
    return cond
end

function affect(integrator)
    terminate!(integrator)
end