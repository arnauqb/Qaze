export save_xi_plots, plot_grid_boundaries, plot_delta_taus
using PyPlot
using Printf
#cmocean = pyimport("cmocean")
LogNorm = matplotlib.colors.LogNorm
function save_xi_plots(xi_grids, wind::WindStruct)
    for i in 1:length(xi_grids)
        fig, ax = plt.subplots()
        ax.set_xlabel("R [Rg]")
        ax.set_ylabel("z [Rg]")
        ax.set_xlim(0, wind.grids.r_range[end])
        ax.set_ylim(0, wind.grids.z_range[end])
        xi_grid = xi_grids[i]
        cm = ax.pcolormesh(wind.grids.r_range, wind.grids.z_range, xi_grid', norm=LogNorm(1e-5, 1e10))#, cmap = cmocean.cm.ice)
        fig.colorbar(cm, ax=ax)
        fig.savefig(@sprintf("animations/xi_grid_%.02d.png", i), dpi=300)
        plt.close("all")
    end
end

function plot_grid_boundaries(wind::WindStruct, ax=nothing, alpha =1.0)
    for r in wind.grids.r_range
        if ax === nothing
            plt.axvline(r, alpha = alpha, color="black")
        else
            ax.axvline(r, alpha = alpha, color="black")
        end
    end
    for z in wind.grids.z_range
        if ax === nothing
            plt.axhline(z, alpha = alpha, color="black")
        else
            ax.axhline(z, alpha = alpha, color="black")
        end
    end
    return nothing
end

function plot_delta_taus(wind::WindStruct, xmin=0, xmax=1000, ymin=0, ymax=1000)
    fig, ax = plt.subplots()
    plot_grid_boundaries(wind, ax, 0.5)
    delta_taus = compute_optical_thickness_grid(wind)
    r_range = copy(wind.grids.r_range)
    z_range = copy(wind.grids.z_range)
    cm = ax.pcolormesh(delta_taus', norm=LogNorm(1e-5, 1))
    #fig.colorbar(cm, ax=ax)
    #ax.set_xlim(xmin, xmax)
    #ax.set_ylim(ymin, ymax)
    #return fig, ax
end