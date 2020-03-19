using PyPlot
using RegionTrees
LogNorm = matplotlib.colors.LogNorm

export plot_cell,
       plot_grid,
       plot_density

function plot_cell(ax, cell)
    v = hcat(collect(vertices(cell.boundary))...)
    x = v[1, [1,2,4,3,1]]
    y = v[2, [1,2,4,3,1]]
    ax.plot(x,y, color="black", alpha=0.3)
end

function plot_grid(quadtree, depth, ax=nothing)
    if ax === nothing
        fig, ax = plt.subplots()
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
    plot_cell.(Ref(ax), children_list)
    return ax
end

function plot_density(wind; xl=0, xh=3000, yl=0, yh=3000, grid=true, depth = 4, nr=1000, nz=1001, vmin=nothing, vmax=nothing)
    fig, ax = plt.subplots()
    grid_den = zeros(Float64, nr, nz)
    r_range = range(xl, stop=xh, length=nr+1)
    z_range = range(yl, stop=yh, length=nz+1)
    for (i, r) in enumerate(r_range[1:end-1])
        for (j, z) in enumerate(z_range[1:end-1])
            leaf = findleaf(wind.quadtree, [r,z])
            line_id = leaf.data[1]
            if line_id == 0
                dens = 1e2
            else
                dens = interpolate_density(wind.lines[line_id], [r,z])
            end
            #dens = leaf.data[1]
            grid_den[i,j] = dens
        end
    end
    cm = ax.pcolormesh(r_range, z_range, grid_den', norm=LogNorm(vmin=vmin, vmax=vmax))
    if grid 
        plot_grid(wind.quadtree, depth, ax)
    end
    fig.colorbar(cm, ax=ax)
    ax.set_xlim(xl,xh)
    ax.set_ylim(yl,yh)
    ax.set_ylabel("z [Rg]")
    ax.set_xlabel("r [Rg]")
    ax.set_title("Density grid")
    return fig, ax
end