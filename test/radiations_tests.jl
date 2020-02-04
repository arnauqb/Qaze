using Qaze
using Test

@testset "nt_rel_factors" begin
    @test nt_rel_factors.([6, 100, 1e8], 0, 6.) ≈ [0, 0.6522661452445981,  0.9960299145167933] rtol=5e-3 atol=0
    @test nt_rel_factors.(10, 0.99, 6.) ≈ 0.22827160704457788 rtol=5e-3 atol=0
    @test nt_rel_factors.(10, 0.99, 10) ≈ 0. rtol=5e-3 atol=0
    @test nt_rel_factors.(500, 0.5, 10) ≈  0.8296494067756125 rtol=5e-3 atol=0
end

@testset "opacity_x" begin
    @test opacity_x(100) == 100
    @test opacity_x(1e6) == 1
end

@testset "X-ray optical depth" begin
    r_range = [10, 10, 50, 100, 900]
    z_range = [5, 5, 20, 800, 5]
    tau_truth = wind.config["wind"]["n_shielding"] * sqrt.(r_range .^2 + z_range .^2) * SIGMA_T * wind.bh.R_g
    f(r,z) = compute_tau_x(r, z, wind)
    @test f.(r_range, z_range) ≈ tau_truth rtol=1e-4 atol=0
end