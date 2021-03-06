using Qaze
using Test

@testset "NT rel factors" begin
    @test nt_rel_factors.([6, 100, 1e8], 0, 6.) ≈ [0, 0.6522661452445981,  0.9960299145167933] rtol=5e-3 atol=0
    @test nt_rel_factors.(10, 0.99, 6.) ≈ 0.22827160704457788 rtol=5e-3 atol=0
    @test nt_rel_factors.(10, 0.99, 10) ≈ 0. rtol=5e-3 atol=0
    @test nt_rel_factors.(500, 0.5, 10) ≈  0.8296494067756125 rtol=5e-3 atol=0
end

@testset "Xray Opacity" begin
    @test opacity_x(100) == 100
    @test opacity_x(1e6) == 1
end

@testset "Ionization Parameter" begin
    r_range = [5, 100, 500, 1000]
    z_range = [10, 50, 1000, 5]
    density = [1e2, 1e4, 1e8, 1e10]
    f(r,z,rho) = ionization_parameter(r, z, rho, 0., wind)
    xi_truth = [332926679880832.3,  33292667988.083233, 33292.66798808323, 416.14794615238657]
    @test f.(r_range, z_range, density) ≈ xi_truth
    @test ionization_parameter(100, 50, 1e8, 1000, wind) <= 1e-20
end

@testset "Force multiplier" begin
    @test force_multiplier(1e5, 1e5, wind) < 1e-3
    @test force_multiplier(1e5, 1e-5, wind) < 1e-3
    @test force_multiplier(1e-8, 1e-5, wind) > 2000
    @test force_multiplier(1e-8, 1e-5, wind) < 3000
    @test force_multiplier(1e-10, 1e-5, wind) ≈ 2300 atol=0 rtol = 0.1
    @test force_multiplier(1e-2, 1e-5, wind) ≈ 6 atol=0 rtol = 0.1
end

@testset "Force radiation" begin
    @test all(isapprox.(integrate(100., 50., wind, include_tau_uv=false, maxevals=10000), [2.4e-6, 2e-6], atol=0, rtol=0.1))
    @test all(isapprox.(integrate(5., 10., wind, include_tau_uv=false, maxevals=10000), [-3.3e-6, 3e-5], atol=0, rtol=0.1))
    #@test all(isapprox.(force_radiation(1000., 0.1, 0., wind, include_tau_uv=false), [1.9e-11, 9e-11], atol=0, rtol=0.2))
end

