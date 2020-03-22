using RegionTrees
import PyPlot
using Printf
LogNorm = PyPlot.matplotlib.colors.LogNorm

export plot_cell,
       plot_grid,
       plot_density,
       plot_linewidth

function plot_cell(ax, cell, xl, xh, yl, yh)
    #if (cell.boundary.origin[1] > xh) || (cell.boundary.origin[2] > yl)
    #    return nothing
    #end
    #if (cell.boundary.origin[1] + cell.boundary.widths[1] < xl) || 
    #    (cell.boundary.origin[2] + cell.boundary.widths[2] < yl)
    #    return nothing
    #end
    v = hcat(collect(vertices(cell.boundary))...)
    x = v[1, [1,2,4,3,1]]
    y = v[2, [1,2,4,3,1]]
    ax.plot(x,y, color="black", alpha=0.3)
end

function plot_grid(quadtree, depth, ax=nothing, xl=0, xh=1000, yl=0, yh=1000)
    if ax === nothing
        fig, ax = PyPlot.plt.subplots()
    end
    children_list = Any[]
    currentchildren = Any[]
    for kid in quadtree.children
        push!(currentchildren, kid)
        push!(children_list, kid)
    end
    currentchildren_aux = Any[]
    dp = depth 
    while dp > 0
        for kid in currentchildren
            if kid.children === nothing
                continue
            end
            for kid2 in kid.children
                push!(children_list, kid2)
                push!(currentchildren_aux, kid2)
            end
        end
        currentchildren = copy(currentchildren_aux)
        currentchildren_aux = []
        dp -= 1
    end
    plot_cell.(Ref(ax), children_list, xl, xh, yl, yh)
    return ax
end

function plot_density(wind; xl=0, xh=3000, yl=0, yh=3000, grid=true, depth = 4, nr=250, nz=251, vmin=nothing, vmax=nothing, ax = nothing, fig=nothing)
    first = false
    if ax === nothing
        first = true
        fig, ax = PyPlot.plt.subplots()
    end
    grid_den = zeros(Float64, nr, nz)
    r_range = range(xl, stop=xh, length=nr+1)
    z_range = range(yl, stop=yh, length=nz+1)
    for (i, r) in enumerate(r_range[1:end-1])
        for (j, z) in enumerate(z_range[1:end-1])
            leaf = findleaf(wind.quadtree, [r,z])
            line_id = leaf.data.line_id
            if line_id == 0
                dens = 1e2#1e2
            else
                dens = interpolate_density(wind.lines[line_id], [r,z], leaf, wind)
            end
            #dens = leaf.data[1]
            grid_den[i,j] = dens
        end
    end
    cm = ax.pcolormesh(r_range, z_range, grid_den', norm=LogNorm(vmin=vmin, vmax=vmax))
    #heatmap(grid_den')
    if grid 
        plot_grid(wind.quadtree, depth, ax, xl, xh, yl, yh)
    end
    if first
        fig.colorbar(cm, ax=ax)
    end
    ax.set_xlim(xl,xh)
    ax.set_ylim(yl,yh)
    ax.set_ylabel("z [Rg]")
    ax.set_xlabel("r [Rg]")
    ax.set_title("Density grid")
    return fig, ax
end

function plot_linewidth(line, ax=nothing)
    if ax === nothing
        fig, ax = PyPlot.plt.subplots()
    end
    r_hist = line.p.u_hist[:,1]
    z_hist = line.p.u_hist[:,2]
    lwhist = line.p.line_width / line.p.r_0 .* r_hist ./ 2
    ax.plot(r_hist, z_hist, "o-", markersize=3, color = "blue")
    ax.plot(r_hist .- lwhist, z_hist, "o-", markersize=3, color = "blue")
    ax.plot(r_hist .+ lwhist, z_hist, "o-", markersize=3, color = "blue")
    return ax
end

function plot_line_gif(line, wind ; xl, xh, yl, yh, depth=6, folder="lineimg", i=0, fig=nothing, ax=nothing)
    #quadtree_initialize(wind)
    k = 1
    while line.sol.retcode == :Default
    #for i in 2:length(line.p.u_hist[:,1])-1 
        step!(line)
        #rmin = minimum(line.p.u_hist[:,1])
        #xl = rmin - rmin * line.p.line_width / line.p.r_0 / 2
        #rmax = maximum(line.p.u_hist[:,1])
        #xh = rmax + rmax * line.p.line_width / line.p.r_0 / 2
        #yl = minimum(line.p.u_hist[:,2])
        #yh = maximum(line.p.u_hist[:,2]) * 2
        lw_norm = line.p.line_width / line.p.r_0
        #quadtree_fill_timestep(currentpoint, previouspoint, n, lw_norm, line.p.line_id, line, line.p.wind)
        fig, ax = plot_density(wind, xl=xl, xh=xh, yl=yl, yh=yh, vmin=nothing, vmax=nothing, grid=false, depth=depth, fig=fig, ax=ax)
        ax.plot(line.p.u_hist[1:k,1], line.p.u_hist[1:k,2], "o-", markersize=3, color="black")
        fig.savefig(folder*@sprintf("/image_%03d.png",i))
        PyPlot.plt.close("all")
        k+=1
        i+=1
    end
    return i, fig, ax
    #end
end

function plot_multiple_lines(i1, i2, wind, folder="gifs", fig=nothing, ax=nothing)
    j = 0
    for i in i1:i2
        wind.lines[i] = initialize_line!(i, wind.lines_range[i], wind)
        rmin = 6#wind.lines[i].p.r_0 - wind.lines[i].p.line_width * 2
        rmax = wind.lines[i].p.r_0 + wind.lines[i].p.line_width * 2
        ymax = 0.1#max(2 * maximum(wind.lines[i].p.u_hist[:,2]), 0.01)
        #folderi = folder*@sprintf("/line_%02d", i)
        j, fig, ax = plot_line_gif(wind.lines[i], wind, folder=folder, xl=rmin, xh=rmax, yl=0, yh=ymax, i=j, fig=fig, ax=ax)
    end
end
