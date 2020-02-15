@testset "Tau UV disk blob" begin
    r_range = range(1., stop=1000., length=20)
    z_range = range(1e-4, stop=1000., length=20)
    r_d_range = range(6., stop=1000., length=20)
    phi_d_range = range(0., stop=pi, length=20)
    for r in r_range
        for z in z_range
            for r_d in r_d_range
                for phi in phi_d_range
                   # println("$r_d, $phi, $r, $z")
                    delta = sqrt(r^2 + z^2 + r_d^2 - 2 * r * r_d * cos(phi))
                    tau = delta * SIGMA_T * wind.grids.n_vacuum * wind.bh.R_g
                    #@test tau_uv_disk_blob(wind, r_d, phi, r, z) ≈ tau atol=0 rtol=1e-2
                    @test tau_uv_disk_blob(r_d, phi, r, z, wind) ≈ tau atol=0 rtol=1e-2
                end
            end
        end
    end
    #deltas = sqrt.(r_range.^2 + r_d_range.^2 + z_range.^2 - 2 .* r_range .* r_d_range .* cos.(phi_d_range))
    #taus = deltas .* SIGMA_T .* wind.grids.n_vacuum .* wind.bh.R_g
    #f(r_d, phi_d, r, z) = tau_uv_disk_blob(wind, r_d, phi_d, r, z)
    #@test all(isapprox.(f.(r_d_range, phi_d_range, r_range, z_range), taus, atol=0, rtol=1e-3))

    #wind.grids.density .*= 4.
    #taus2 = 4 .* taus
    #@test all(isapprox.(f.(r_d_range, phi_d_range, r_range, z_range), taus2, atol=0, rtol=1e-3))
    #wind.grids.density ./= 4.

    #wind.grids.fm .= 10
    #taus3 = (1 + 10) .* taus
    #@test all(isapprox.(f.(r_d_range, phi_d_range, r_range, z_range), taus3, atol=0, rtol=1e-3))
    #wind.grids.fm .= 0
end

@testset "integrands" begin
    v_base = [0., 0.]
    v = [0., 0.]
    integrate_notau_kernel(v_base, 100, 1, 200, 50, wind)
    wind.grids.mdot .*= 2
    integrate_notau_kernel(v, 100, 1, 200, 50, wind)
    @test all(isapprox.(v, 2*v_base, atol=0, rtol=1e-4))
    wind.grids.mdot ./= 2
    wind.grids.uv_fractions .*= 2
    integrate_notau_kernel(v, 100, 1, 200, 50, wind) 
    @test all(isapprox.(v, 2*v_base, atol=0, rtol=1e-4))
    wind.grids.uv_fractions ./= 2

    integrate_kernel(v_base, 100, 1, 200, 50, wind)
    wind.grids.mdot .*= 2
    integrate_kernel(v, 100, 1, 200, 50, wind)
    @test all(isapprox.(v, 2*v_base, atol=0, rtol=1e-4))
    wind.grids.mdot ./= 2
    wind.grids.uv_fractions .*= 2
    integrate_kernel(v, 100, 1, 200, 50, wind) 
    @test all(isapprox.(v, 2*v_base, atol=0, rtol=1e-4))
    wind.grids.uv_fractions ./= 2
    #wind.grids.density .= 1e20
    #integrate_kernel(v, 100, 1, 200, 50, wind) 
    #@test all(v .< [1e-11,1e-11])
    #wind.grids.density .= wind.grids.n_vacuum
    #wind.grids.density .= 1e10
    #integrate_kernel(v, 100, 1, 200, 50, wind) 
    #@test all(v .< [1e-11,1e-11])
    #wind.grids.density .= wind.grids.n_vacuum
end