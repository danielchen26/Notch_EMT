using Makie
using GLMakie
Base.@kwdef mutable struct Lorenz
    dt::Float64 = 0.01
    σ::Float64 = 10
    ρ::Float64 = 28
    β::Float64 = 8 / 3
    x::Float64 = 1
    y::Float64 = 1
    z::Float64 = 1
end

function step!(l::Lorenz)
    dx = l.σ * (l.y - l.x)
    dy = l.x * (l.ρ - l.z) - l.y
    dz = l.x * l.y - l.β * l.z
    l.x += l.dt * dx
    l.y += l.dt * dy
    l.z += l.dt * dz
    Point3f(l.x, l.y, l.z)
end

attractor = Lorenz()

points = Node(Point3f[])
colors = Node(Int[])

set_theme!(backgroundcolor = :white)

fig, ax, l = lines(points, color = colors,
    colormap = :inferno, transparency = true,
    axis = (; type = Axis3, protrusions = (0, 0, 0, 0),
        viewmode = :fit, limits = (-30, 30, -30, 30, 0, 50)))

record(fig, "lorenz.mp4", 1:120) do frame
    for i = 1:50
        push!(points[], step!(attractor))
        push!(colors[], frame)
    end
    ax.azimuth[] = 1.7pi + 0.3 * sin(2pi * frame / 120)
    notify.((points, colors))
    l.colorrange = (0, frame)
end




## =============
x = range(0, 10, length = 100)
y = sin.(x)
lines(x, y)

using CairoMakie
CairoMakie.scatter(x, y)


x = range(0, 10, length = 100)
y1 = sin.(x)
y2 = cos.(x)

lines(x, y1)
lines!(x, y2)
current_figure()

CairoMakie.scatter(x, y1, markersize = 5, colormap = :thermal, colorrange = (-1, 1), label = "sin")
axislegend()
current_figure()




fig = Figure(resolution = (2560, 1440))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
ax3 = Axis(fig[2, 1:2])
fig

f(x) = x^3 - 2x^2 + 3x - 1
fig = Figure()
ax1, l1 = lines(fig[1, 1], 0 .. 10, f, color = :red)
ax2, l2 = lines(fig[2, 1], 0 .. 10, cos, color = :blue)
Legend(fig[1:2, 2], [l1, l2], ["f(x)", "cos"])
fig
fig = Figure(resolution = (800, 600))
ax, hm = CairoMakie.heatmap(fig, randn(20, 20))
Colorbar(fig[1, 2], hm)
fig