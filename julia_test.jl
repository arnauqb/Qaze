using PyPlot
LogNorm = matplotlib.colors.LogNorm
function plot()
    fig, ax = plt.subplots()
    x_range = range(0, stop=100, length = 100)
    y_range = range(0, stop=100, length = 101)
    grid = rand(100, 101)
    cm = ax.pcolormesh(x_range, y_range, grid', norm=LogNorm())
    return fig, ax
end