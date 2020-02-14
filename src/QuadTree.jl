using RegionTrees: vertices, Cell, findleaf, children, split!
import Base: copy, isequal, ==
using StaticArrays: SVector
using RegionTrees
export compute_cell_intersection, compute_cell_size, fill_linewidth, fill_point, 
fill_all_lines, copy, refine_all, refine_leaf, compute_density_grid, 
plot_density_grid_tree, plot_tricontour, plot_taux_grid_tree, fill_line, 
fill_and_refine, fill_and_refine_linewidth, fill_and_refine_line, fill_and_refine_all_lines, erase_line_from_tree
using PyPlot
LogNorm = matplotlib.colors.LogNorm

function compute_cell_size(cell::Cell)
    r1, z1 = vertices(cell.boundary)[1,1]
    r2, z2 = vertices(cell.boundary)[2,2]
    distance = sqrt((r2-r1)^2 + (z2-z1)^2)
    return distance
end

function copy(cell::Cell)
    newcell = Cell(cell.boundary, cell.data)
    return newcell
end

function isequal(cell1::Cell, cell2::Cell)
    return cell1.boundary == cell2.boundary
end

function ==(cell1::Cell, cell2::Cell)
    return cell1.boundary == cell2.boundary
end

function reinitialize_tree(wind::WindStruct)
    quadtree_max_radius = wind.config["grids"]["r_max"]
    quadtree_max_height = wind.config["grids"]["z_max"]
    n_vacuum = wind.grids.n_vacuum
    delta_tau_0 = n_vacuum * sqrt(quadtree_max_height^2 + quadtree_max_radius^2) * wind.bh.R_g * SIGMA_T
    wind.quadtree = Cell(SVector(0., 0.), SVector(2 * quadtree_max_radius, 2* quadtree_max_height), [n_vacuum, 0., delta_tau_0])
end

function compute_cell_intersection(current_point, cell::Cell, initial_point, final_point)
    r, z = current_point
    r_i, z_i = initial_point
    r_f, z_f = final_point
    direction_r = sign(r_f - r_i)
    direction_z = sign(z_f - z_i)
    cell_vertices = vertices(cell)
    if direction_r == 1
        cell_r_boundary = cell_vertices[2,2][1]
    else
        cell_r_boundary = cell_vertices[1,1][1]
    end
    if direction_z == 1
        cell_z_boundary = cell_vertices[2,2][2]
    else
        cell_z_boundary = cell_vertices[1,1][2]
    end
    if r_i == r_f
        lambda_r = Inf
    else
        lambda_r = (cell_r_boundary - r) / (r_f - r_i)
    end
    if z_i == z_f
        lambda_z = Inf
    else
        lambda_z = (cell_z_boundary - z) / (z_f - z_i)
    end
    lambda = min(lambda_r, lambda_z)
    try
        @assert lambda >= 0.0
    catch
        println("current point : $current_point")
        println("cell: $cell")
        println("initial point: $initial_point")
        println("final point: $final_point")
        throw(DomainError)
    end
    intersection = [r, z] + lambda .* [r_f - r_i, z_f - z_i]
    return intersection
end

function fill_linewidth(point, linewidth, density, fm, wind::WindStruct)
    r_min = wind.grids.r_range[1]
    r_max = wind.grids.r_range[end]
    point1 = [max(point[1] - linewidth / 2., r_min) , point[2]] 
    point2 = [min(point[1] + linewidth / 2., r_max) , point[2]] 
    sigmarg = SIGMA_T * wind.bh.R_g
    currentpoint = copy(point1)
    currentleaf = findleaf(wind.quadtree, point1)
    previousleaf = copy(currentleaf)
    point2leaf = findleaf(wind.quadtree, point2)
    n_fill = max(density, currentleaf.data[1])
    fm_fill = max(fm, currentleaf.data[2])
    deltad = compute_cell_size(currentleaf)
    deltatau = sigmarg * deltad * n_fill
    currentleaf.data = [n_fill, fm_fill, deltatau]
    while (currentleaf != point2leaf)
        currentpoint = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        if currentleaf.boundary == previousleaf.boundary
            currentpoint .-= 1e-8
            currentleaf = findleaf(wind.quadtree, currentpoint)
        end
        n_fill = max(density, currentleaf.data[1])
        fm_fill = max(fm, currentleaf.data[2])
        deltad = compute_cell_size(currentleaf)
        deltatau = sigmarg * deltad * n_fill
        currentleaf.data = [n_fill, fm_fill, deltatau]
        previousleaf = copy(currentleaf)
    end
end

function fill_point(point1, point2, linewidth_normalized, density, fm, wind::WindStruct)
    if point1[1] < point2[1]
        point1, point2 = point2, point1
    end
    currentleaf = findleaf(wind.quadtree, point1)
    previousleaf = copy(currentleaf)
    point2leaf = findleaf(wind.quadtree, point2)
    linewidth = linewidth_normalized * point1[1]
    fill_linewidth(point1, linewidth, density, fm, wind)
    currentpoint = copy(point1)
    while currentleaf != point2leaf
        #println("point1 : $currentpoint, point2 : $point2, currentleaf : $currentleaf")
        currentpoint = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        if currentleaf.boundary == previousleaf.boundary
            #println("currentleaf: $(currentleaf.boundary)")
            currentpoint .-= 1e-8
            #println(currentpoint)
            currentleaf = findleaf(wind.quadtree, currentpoint)
            #println("newleaf: $(currentleaf.boundary)")
        end
        linewidth = linewidth_normalized * currentpoint[1]
        fill_linewidth(currentpoint, linewidth, density, fm, wind)
        previousleaf = copy(currentleaf)
    end
end

function fill_line(line, wind::WindStruct)
    if length(line.p.n_hist) < 3
        return nothing
    end
    linewidth_normalized = line.p.line_width / line.p.r_0 
    for i in 2:(length(line.p.n_hist) - 1)
        point1 = line.p.u_hist[i, 1:2]
        point2 = line.p.u_hist[i+1, 1:2]
        rho = line.p.n_hist[i]
        fm = line.p.fm_hist[i]
        fill_point(point1, point2, linewidth_normalized, rho, fm, wind)
    end
end

function fill_all_lines(wind::WindStruct)
    #initialize them all
    for leaf in allleaves(wind.quadtree)
        leaf.data = [wind.grids.n_vacuum, 0., 0.01]
    end
    k=1
    for line in wind.lines
        #println("filling line ... $k")
        k += 1
        fill_line(line, wind)
    end
end

function split_cell_initialization(cell, child_indices)
    newdata = [1e2, 0., 0.] 
    return newdata
end

function refine_leaf(currentpoint, leaf, n, fm, wind)
    optically_thin = false
    while !optically_thin
        optically_thin = true
        maxz = vertices(leaf.boundary)[2,2][2]
        refine_condition = (leaf.data[3] > 0.1) && (maxz > wind.config["wind"]["z_0"]) 
        if refine_condition
            optically_thin = false
            split!(leaf, split_cell_initialization)
            newleaf = findleaf(leaf, currentpoint)
            deltad = compute_cell_size(newleaf)
            deltatau = n * SIGMA_T * wind.bh.R_g * deltad
            newleaf.data = [n, fm, deltatau]
            leaf = newleaf
        end
    end
end

function fill_and_refine_linewidth(point, linewidth, density, fm, wind::WindStruct)
    r_min = 0 
    r_max = 6000 
    point1 = [max(point[1] - linewidth / 2., r_min) , point[2]] 
    point2 = [min(point[1] + linewidth / 2., r_max) , point[2]] 
    sigmarg = SIGMA_T * wind.bh.R_g
    currentpoint = copy(point1)
    currentleaf = findleaf(wind.quadtree, point1)
    point2leaf = findleaf(wind.quadtree, point2)
    n_fill = max(density, currentleaf.data[1])
    fm_fill = max(fm, currentleaf.data[2])
    deltad = compute_cell_size(currentleaf)
    deltatau = sigmarg * deltad * n_fill
    currentleaf.data = [n_fill, fm_fill, deltatau]
    refine_leaf(currentpoint, currentleaf, n_fill, fm_fill, wind)
    currentleaf = findleaf(wind.quadtree, currentpoint)
    point2leaf = findleaf(wind.quadtree, point2)
    while (currentleaf != point2leaf)
        currentpoint = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        n_fill = max(density, currentleaf.data[1])
        fm_fill = max(fm, currentleaf.data[2])
        deltad = compute_cell_size(currentleaf)
        deltatau = sigmarg * deltad * n_fill
        currentleaf.data = [n_fill, fm_fill, deltatau]
        refine_leaf(currentpoint, currentleaf, n_fill, fm_fill, wind)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        point2leaf = findleaf(wind.quadtree, point2)
    end
end

function fill_and_refine(point1, point2, linewidth_normalized, n, fm, wind::WindStruct)
    lw = linewidth_normalized * point1[1]
    point1leaf = findleaf(wind.quadtree, point1)
    currentpoint = copy(point1)
    currentleaf = copy(point1leaf)
    previousleaf = copy(currentleaf)
    point2leaf = findleaf(wind.quadtree, point2)
    fill_and_refine_linewidth(currentpoint, lw, n, fm, wind)
    currentleaf = findleaf(wind.quadtree, currentpoint)
    point2leaf = findleaf(wind.quadtree, point2)
    while(currentleaf != point2leaf)
        try
            currentpoint = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        catch
            println("currentpoint : $currentpoint")
            println("currentleaf: $currentleaf")
            println("point1: $point1")
            println("point2: $point2")
            println("point1leaf: $point1leaf")
            println("point2leaf: $point2leaf")
            throw(DomainError)
        end
        currentleaf = findleaf(wind.quadtree, currentpoint)
        if currentleaf == previousleaf
            currentpoint .-= 1e-8
            currentleaf = findleaf(wind.quadtree, currentpoint)
        end
        fill_and_refine_linewidth(currentpoint, lw, n, fm , wind)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        point2leaf = findleaf(wind.quadtree, point2)
        previousleaf = copy(currentleaf) 
    end
end

function fill_and_refine_line(line, wind::WindStruct)
    for i in 1:length(line.p.n_hist)-1
        point1 = line.p.u_hist[i, 1:2]
        point2 = line.p.u_hist[i+1, 1:2]
        linewidth_normalized = line.p.line_width / line.p.r_0
        n = line.p.n_hist[i]
        fm = line.p.fm_hist[i]
        fill_and_refine(point1, point2, linewidth_normalized, n, fm, wind)
    end
end

function fill_and_refine_all_lines(wind::WindStruct)
    for i in 1:length(wind.lines)
        if !isassigned(wind.lines, i)
            continue
        end
        line = wind.lines[i]
        fill_and_refine_line(line)
    end
end

function refine_all(wind::WindStruct)
    fill_all_lines(wind)
    optically_thin = false
    while (!optically_thin)
        println("refining ...")
        optically_thin = true
        for leaf in allleaves(wind.quadtree)
            maxz = vertices(leaf.boundary)[2,2][2]
            refine_condition = (leaf.data[3] > 0.1) && (maxz > wind.config["wind"]["z_0"]) 
            if refine_condition
                optically_thin = false
                split!(leaf)
            end
        end
        fill_all_lines(wind)
    end
end

function compute_density_grid(wind::WindStruct)
    grid = zeros(Float64, 2000, 2001)
    r_range = range(wind.grids.r_range[1], stop=wind.grids.r_range[end], length=2000)
    z_range = range(wind.grids.z_range[1], stop=wind.grids.z_range[end], length=2001)
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            leaf = findleaf(wind.quadtree, [r,z])
            dens = leaf.data[1]
            grid[i,j] = dens
        end
    end
    return r_range, z_range, grid
end

function plot_taux_grid_tree(wind::WindStruct)
    plt.figure()
    nr = 50
    nz = 51
    grid = zeros(Float64, nr, nz)
    r_range = range(wind.grids.r_range[1], stop=wind.grids.r_range[end], length=nr)
    z_range = range(wind.config["wind"]["z_0"], stop=wind.grids.z_range[end], length=nz)
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            taux = compute_taux_tree(r, z, wind)
            grid[i,j] = taux
        end
    end
    plt.pcolormesh(r_range, z_range, grid', norm=LogNorm(vmax=1e1))
    plt.colorbar()
    display(gcf())
end

function plot_density_grid_tree(wind::WindStruct)
    plt.figure()
    nr = 500
    nz = 501
    grid = zeros(Float64, nr, nz)
    r_range = range(wind.grids.r_range[1], stop=wind.grids.r_range[end], length=nr)
    z_range = range(wind.config["wind"]["z_0"], stop=wind.grids.z_range[end], length=nz)
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            leaf = findleaf(wind.quadtree, [r,z])
            dens = leaf.data[1]
            grid[i,j] = dens
        end
    end
    plt.pcolormesh(r_range, z_range, grid', norm=LogNorm(vmin=1e5, vmax=1e8))
    plt.colorbar()
    println("printing grid...")
    for leaf in allleaves(wind.quadtree)
        v = hcat(collect(vertices(leaf.boundary))...)
        x = v[1, [1,2,4,3,1]]
        y = v[2, [1,2,4,3,1]]
        plt.plot(x,y, color="black", alpha=0.3)
    end
    #for line in wind.lines
    #    plt.plot(line.p.u_hist[:,1], line.p.u_hist[:,2], color="white", linewidth=0.5)
    #end
    #display(gcf())
end

function plot_tricontour(wind::WindStruct)
    r_range = []
    z_range = []
    rho_range = []
    for leaf in allleaves(wind.quadtree)
        r, z = center(leaf)
        rho = leaf.data[1]
        push!(r_range, r)
        push!(z_range, z)
        push!(rho_range, rho)
    end
    plt.figure()
    #plt.tripcolor(r_range, z_range, rho_range, norm=LogNorm())
    plt.tricontourf(r_range, z_range, rho_range, 20, norm=LogNorm())
end

function erase_line_from_tree(line_id, wind::WindStruct)
    println("resetting tree...")
    reinitialize_tree(wind)
    line = wind.lines[line_id]
    n_vacuum = wind.grids.n_vacuum
    line.n_hist .= n_vacuum * ones(Float64, size(line.n_hist))
    line.fm_hist .= zeros(Float64, size(line.fm_hist))
    println("refilling and refining all lines...")
    fill_and_refine_all_lines(wind)
end