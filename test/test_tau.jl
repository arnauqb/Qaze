using Test
using ProgressMeter
using RegionTrees
export density_profile,
       density_profile_simple,
       tau_densityprofile_analytical,
       fill_constant_density_grid,
       fill_densityprofile_grid,
       fill_densityprofile_simple_grid,
       tau_densityprofile_analytical,
       tau_densityprofile_analytical_simple

wind = initialize(string(@__DIR__,"/config_tau_tests.toml"))

function density_profile(r, z, n0, r0, z0)
    return n0 * r / r0 * exp(-z / z0)
end

function density_profile_simple(z, n0, r0, z0)
    return n0 * z / z0
end

function fill_constant_density_grid(wind, n0)
    quadtree_initialize(wind)
    split!(wind.quadtree, dumb_fill_data)
    needs_refine = true
    while needs_refine 
        needs_refine = false
        for leaf in allleaves(wind.quadtree)
            if cell_width(leaf) > 1
                needs_refine = true
                split!(leaf, dumb_fill_data)
            end
        end
    end
    for leaf in allleaves(wind.quadtree)
        r_range = 10 .^ range(log10(leaf.boundary.origin[1] + 1e-6), 
                              log10(leaf.boundary.origin[1] + leaf.boundary.widths[1]), length=3)[2:end-1]
        z_range = 10 .^ range(log10(leaf.boundary.origin[2] + 1e-6), 
                              log10(leaf.boundary.origin[2] + leaf.boundary.widths[2]), length=10)[2:end-1]
        for r in r_range
            for z in z_range
                density = n0 
                quadtree_fill_point([r,z], 1, density, 1, nothing, wind, needs_refinement=false)
            end
        end
    end
end

function dumb_fill_data(cell, child_indices)
    return CellData(Int[], Float64[], Float64[], Float64[])
end

function fill_densityprofile_simple_grid(wind, n0, r0, z0)
    quadtree_initialize(wind)
    split!(wind.quadtree, dumb_fill_data)
    needs_refine = true
    println("Refining...")
    while needs_refine 
        needs_refine = false
        for leaf in allleaves(wind.quadtree)
            if cell_width(leaf) > 0.2
                needs_refine = true
                split!(leaf, dumb_fill_data)
            end
        end
    end
    println("Done")
    println("Filling...")
    for leaf in allleaves(wind.quadtree)
        r_range = 10 .^ range(log10(leaf.boundary.origin[1] + 1e-6), 
                              log10(leaf.boundary.origin[1] + leaf.boundary.widths[1]), length=3)[2:end-1]
        z_range = 10 .^ range(log10(leaf.boundary.origin[2] + 1e-6), 
                              log10(leaf.boundary.origin[2] + leaf.boundary.widths[2]), length=10)[2:end-1]
        for r in r_range
            for z in z_range
                density = density_profile_simple(z, n0, r0, z0)
                quadtree_fill_point([r,z], 1, density, 1, nothing, wind, needs_refinement=false)
            end
        end
    end
    #r_range = range(0, 100, length=500)
    #z_range = 10 .^ range(-6, 2, length=1000)
    #for z in z_range
    #    density = density_profile_simple(z, n0, r0, z0)
    #    #quadtree_fill_point([r,z], 0.1, density, 1, nothing, wind)
    #    quadtree_fill_horizontal(0, 100, z, 0.1, density, 1, nothing, wind, true)
    #end
end

function tau_densityprofile_analytical(rd, r, z, n0, r0, z0)
    term1 = sqrt(((r - rd)^2 + z^2)) * n0 / r0 
    term2 = (r-rd) * (z0/z)^2 * (1- (1 + z / z0) * exp(-z/z0))
    term3 = rd * (z0 / z) * (1 - exp(-z/z0))
    return term1 * (term2 + term3)
end

function tau_densityprofile_analytical_simple(rd, r, z, n0, z0)
    term1 = 0.5 * sqrt((r-rd)^2 + z^2)
    term2 = n0 / z0 * z
    return term1 * term2
end

#@testset "Test constant density grid" begin
#    density = 1e8
#    wind.config["grids"]["r_max"] = 256.0
#    wind.config["grids"]["z_max"] = 256.0
#    fill_constant_density_grid(wind, density)
#    sigmarg = wind.bh.R_g * SIGMA_T
#    r_d_range = range(0, 100, length=50)
#    r_range = range(0, 250, length=50)
#    z_range = 10 .^ range(-9, log10(250), length=50)
#    @showprogress 1 "Computing taus..." for rd in r_d_range
#        for r in r_range
#            for z in z_range
#                tautrue = sqrt((r-rd)^2 + z^2) * sigmarg * density
#                tau = tau_uv_disk_blob(rd, r, z, wind, Inf) * sigmarg
#                res = isapprox(tautrue, tau, atol=0, rtol=1e-2)
#                if !res 
#                    println("rd: $rd")
#                    println("r: $r")
#                    println("z: $z")
#                    println("tautrue: $tautrue")
#                    println("tau : $tau")
#                end
#                @test res
#            end
#        end
#    end
#end

@testset "Test simple profile density grid" begin
    n0 = 1e8
    r0 = 10
    z0 = 1
    wind.config["r_max"] = 60.0
    wind.config["z_max"] = 6.0
    fill_densityprofile_simple_grid(wind, n0, r0, z0)
    sigmarg = wind.bh.R_g * SIGMA_T
    r_d_range = range(0, 90, length=10)
    r_range = range(0, 90, length=10)
    z_range = 10 .^ range(-5, 0.5, length=10)
    @showprogress 1 "Computing taus..." for rd in r_d_range
        for r in r_range
            for z in z_range
                tautrue = tau_densityprofile_analytical_simple(rd, r, z, n0, r0, z0) 
                tau = 0.
                try
                    tau = tau_uv_disk_blob(rd, r, z, wind, Inf)
                catch
                    println("rd: $rd")
                    println("r: $r")
                    println("z: $z")
                    throw(DomainError)
                end
                res = isapprox(tautrue, tau, atol=0, rtol=1)
                if !res 
                    println("rd: $rd")
                    println("r: $r")
                    println("z: $z")
                    println("tautrue: $tautrue")
                    println("tau : $tau")
                end
                @test res
            end
        end
    end
end


#@testset "Test profile density grid" begin
#    n0 = 1e8
#    r0 = 10
#    z0 = 1
#    fill_densityprofile_grid(wind, n0, r0, z0)
#    sigmarg = wind.bh.R_g * SIGMA_T
#    r_d_range = range(0, 100, length=50)
#    r_range = range(0, 250, length=50)
#    z_range = 10 .^ range(-9, log10(250), length=50)
#    @showprogress 1 "Computing taus..." for rd in r_d_range
#        for r in r_range
#            for z in z_range
#                density = density_profile(r, z, n0, r0, z0)
#                tautrue = sqrt((r-rd)^2 + z^2) * sigmarg * density
#                tau = tau_uv_disk_blob(rd, r, z, wind, Inf) * sigmarg
#                @test isapprox(tautrue, tau, atol=0, rtol=1e-2)
#            end
#        end
#    end
#end