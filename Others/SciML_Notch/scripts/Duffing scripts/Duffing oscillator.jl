##
using Symbolics
using Plots
using DifferentialEquations
using ProgressMeter, Measures
using DataFrames
using LaTeXStrings
using IntervalArithmetic, IntervalRootFinding
cloud = "/Users/chentianchi/Cubbit/Notch_EMT/"
function add_Boundary_event(switch_amplitude, switch_frequency, plt; lw = 2)
  df = DataFrame(switch_amplitude = switch_amplitude, switch_frequency = switch_frequency)
  # CSV.write("switching_freq_amplitude.csv",df)

  gp = groupby(df, :switch_frequency)
  boundary_amplitude = []
  for i in gp
    append!(boundary_amplitude, minimum(i.switch_amplitude))
  end

  critical_freq = [keys(gp)[i].switch_frequency for i = 1:gp.ngroups]
  boundary_amplitude
  plt_boundary = plot!(plt, boundary_amplitude, critical_freq, label = "switching boundary", lw = lw)
end

## ============ adding functions =========


##  parameters
A = 0.1;
ω = 1.4;
c = -20;
β = 10;
γ = -1;

# function
V(x) = 1/2*c*x^2 + 1/4*β*x^4   # potential function
x_pos = sqrt(-c / β)               # Eq position
E = c^2 / (4 * β)                    # magnitude of energy barrier
V(x_pos)                         # energy relative to zero level

tt = -10:0.1:10
plot(tt, V.(tt), ylim = [-20,50])

## ------------- linear case ---------
γ = .1
c = 10

A = 10.5
ω = .5

β = 0   # linear case β = 0
θ = 0 # linear case θ = 0
ϕ = 1 # signal phase
# ------------ potential function (asymmetric)
V(x) = θ * x + 1 / 2c * x^2 + 1 / 4β * x^4
Vp(x) = c * x + β * x^3 + θ

function duffing(du, u, p, t)
  x, v = u
  c, β, γ, A, ω = p
  du[1] = u[2]
  du[2] = -c * u[1] - β * u[1]^3 - γ * u[2] + A * (1 + sign(cos(ω * t + ϕ))) # pulsatile signal
end


u0 = [0.0; 0.0]
p = [c, β, γ, A, ω, ϕ]
Δt = 1e-2
T = 20
tspan = (0.0, T)
prob = ODEProblem(duffing, u0, tspan, p)
sol = solve(prob, RK4(), adaptive = false, dt = Δt)

rts = values(IntervalRootFinding.roots(x -> sol(x)[1], 0 .. 10))
midpoints = mid.(interval.(rts))
@show sort(midpoints)

# Visualization
plotly()
plt_ds1 = plot(sol,  idxs = [1], label = ["Position" "Velocity"], lw = [0.5 1])
hline!([0], linestyle=:dash)
ylims!(-6,6)
xticks!(0:1:20)
title!("freq is $ω")


# plt_test = plot(plt_ds1,plt_ds2,plt_ds3)


##
# pos = sol[1, :]
# velocity = sol[2, :]
# E_tot = 1 / 2 * velocity .^ 2 .+ V.(pos)
# plt_E = plot(sol.t, E_tot, label = "Energy")
#
# plt = plot(plt_ds, plt_E, layout = (2, 1), dpi = 300)
#
# # ======  potential surface
# xs = -5:0.1:5
# plt_V = plot(xs, V.(xs), title = "Potential Surface", dpi = 300)
#
# l = @layout [a; b c]
# plot(plt_ds, plt_E, plt_V, layout = l, dpi = 300)
##

##
function linear(; c = 10, β = 0, γ = 0.3, A_range = 0.0:50, ω_range = 0:0.2:30)
  Aw_set = []
  # c = 10; β= 0; γ = 0.3;
  c = c
  β = β
  γ = γ
  @showprogress for ω = ω_range#collect(exp10.(-4:0.1:2))
    for A = A_range#collect(exp10.(-1:0.1:3))
      p[1] = c
      p[2] = β
      p[3] = γ
      p[4] = A
      p[5] = ω
      p[6] = ϕ
      # @show A, ω
      prob1 = remake(prob, p = p)
      sol = solve(prob1, RK4(), adaptive = false, dt = Δt)
      plt_ds = plot(sol, label = ["Position" "Velocity"], lw = [3 1])
      # display(plt_ds)
      pos = sol[1, :]
      velocity = sol[2, :]
      # @show maximum(abs.(pos))
      pos_later = pos[length(sol.t), end] # try to get the steady state solution max.
      cri = maximum(abs.(pos_later)) > 1.0
      cri ? push!(Aw_set, [A, ω]) : nothing
      if cri
        break
      end
  
    end
  end
  #Plot A w curve
  MM = hcat(Aw_set...)'
  MA = MM[:, 1]
  Mw = MM[:, 2]
  # scatter(MA, Mw,
  #         xlabel = "Wwitching Amplitude", ylabel = "Switching Frequency",
  #         label = "swithing event",
  #         dpi = 500)


  plt = scatter(MA, Mw,
    xlabel = "Switching Amplitude (A)", ylabel = "Switching Frequency (ω)",
    m = (:viridis, 0.4, Plots.stroke(2, :green)),
    label = "swithing event", dpi = 500)
  add_Boundary_event(MA, Mw, plt, lw = 1.5)
  return plt
end



function linear_ϕ_independent(; c = 10, β = 0, γ = 0.3, A_range = 0.0:50, ω_range = 0:0.2:30)
  Aw_set = []
  # c = 10; β= 0; γ = 0.3;
  c = c
  β = β
  γ = γ
  @showprogress for ω = ω_range#collect(exp10.(-4:0.1:2))
                  for A = A_range#collect(exp10.(-1:0.1:3))
                    switch_set = []
                    for ϕ = 0:2π
                      p[1] = c
                      p[2] = β
                      p[3] = γ
                      p[4] = A
                      p[5] = ω
                      p[6] = ϕ
                      # @show A, ω
                      prob1 = remake(prob, p = p)
                      sol = solve(prob1, RK4(), adaptive = false, dt = Δt)
                      plt_ds = plot(sol, label = ["Position" "Velocity"], lw = [3 1])
                      # display(plt_ds)
                      pos = sol[1, :]
                      velocity = sol[2, :]
                      # @show maximum(abs.(pos))
                      pos_later = pos[length(sol.t), end] # try to get the steady state solution max.
                      cri = maximum(abs.(pos_later)) > 1.0
                      cri ? push!(switch_set, -1) : push!(switch_set, 1)
                      # cri ? push!(Aw_set, [A, ω]) : nothing
                      if cri
                        break
                      end
                    end
                    push!(Aw_set, [A, ω])
                  end
                end
  #Plot A w curve
  MM = hcat(Aw_set...)'
  MA = MM[:, 1]
  Mw = MM[:, 2]
  # scatter(MA, Mw,
  #         xlabel = "Wwitching Amplitude", ylabel = "Switching Frequency",
  #         label = "swithing event",
  #         dpi = 500)


  plt = scatter(MA, Mw,
    xlabel = "Switching Amplitude (A)", ylabel = "Switching Frequency (ω)",
    m = (:viridis, 0.4, Plots.stroke(2, :green)),
    label = "swithing event", dpi = 500)
  add_Boundary_event(MA, Mw, plt, lw = 1.5)
  return plt
end







@showprogress for γ in 0:0.5:1
  plt = linear(; c = 10, β = 0, γ = γ, A_range = 0.0:50, ω_range = 0:0.2:30)
  # plt = linear_ϕ_independent(; c = 10, β = 0, γ = γ, A_range = 0.0:10:50, ω_range = 0:2:10)
  # save_path = joinpath(pwd(), "figures", "linear_pulsatile_steady_lr_ϕ_indep/") # generate path to save
  save_path = joinpath(cloud, "figures", "linear_pulsatile_steady_lr_ϕ_indep/") # generate path to save
  isdir(save_path) || mkpath(save_path)
  savefig(plt, save_path * "c = 10, γ = $γ.png")
end



## 🍏 ================ working on right now ================
γ = 0.1
c = -10
β = 1
A = 5.5
ω = 0.8
ϕ = 0
V(x) = 1/2*c*x^2 + 1/4*β*x^4 
x_pos = sqrt(-c / β)
E_barrier = c^2 / (4 * β) # energy barrier
u0 = [x_pos; 0.0]
p = [c, β, γ, A, ω, ϕ]
Δt = 1e-2
T = 50
tspan = (0.0, T)
# ---- plot potential well --------------------------------
plt_potential = plot(-2x_pos:0.1:2x_pos, V.(-2x_pos:0.1:2x_pos), legend = false)
hline!([-E_barrier],linestyle=:dash)
ylims!(-1.2E_barrier,1.2E_barrier)
xticks!(-2floor(x_pos):1:2floor(x_pos))
title!("Energy barrier is $E_barrier")

function duffing(du, u, p, t)
  x, v = u
  c, β, γ, A, ω, ϕ = p
  du[1] = u[2]
  du[2] = -c * u[1] - β * u[1]^3 - γ * u[2] + A * (1 + sign(cos(ω * t + ϕ))) # pulsatile signal
end

prob = ODEProblem(duffing, u0, tspan, p)
sol = solve(prob, RK4(), adaptive = false, dt = Δt)
plt = plot(sol,label = ["Position" "Velocity"], lw = [3 1])

plot(plt_potential,plt, layout = (2,1))

##
Aw_set = []
suptitle = plot(title = "Duffing Oscillator", framestyle = nothing, showaxis = false, xticks = false, yticks = false, margin = 0mm)
# c = -20;
# β = 1;
# γ = .5;

record_fpt =[]
@showprogress for ω = collect(0:0.05:10)
  for A = 50#collect(exp10.(-1:0.1:3))
    A = round(A, digits =2)
    ω = round(ω, digits =3)
    # Initial totoal energy
    E_init = V(x_pos) # because initial condition v_init = 0
    E_barrier = c^2 / (4 * β) # energy barrier
    p[1] = c
    p[2] = β
    p[3] = γ
    p[4] = A
    p[5] = ω
    p[6] = ϕ
    @show A, ω
    prob1 = remake(prob, p = p)
    sol = solve(prob1, RK4(), adaptive = false, dt = Δt)

    # ---- get zeros point, which means crossing the energy barrier
    rts = values(IntervalRootFinding.roots(x -> sol(x)[1], 0 .. T))
    midpoints = mid.(interval.(rts))
    @show sort(midpoints)
    if isempty(midpoints) == false
      push!(record_fpt,(first(sort(midpoints)), ω))
    end
    # ---- record the first time system 
    ### ==== Visualization ======

    # ===== Dynamics plot
    plt_ds = plot(sol, label = ["Position" "Velocity"],
      xlabel = "Time (t)", color = ["blue" "orange"], linewidth = [2.5 1.5],
      title = "Dynamics (A: $A, ω = $ω)", titlefontsize = 10)
    lb, up = extrema(sol[1, :])
    lim = maximum([abs(lb),abs(up)])
    ylims!(-lim, lim)
    # display(plt_ds)
    # # ===== Total energy plot
    pos = sol[1, :]
    velocity = sol[2, :]
    E_k = 1 / 2 * velocity .^ 2
    E_p = V.(pos)
    E_tot = E_k .+ E_p

    plt_Ek = plot(sol.t, E_k, color = "orange", label = "Kinetic Energy")
    plt_Ep = plot(sol.t, E_p, color = "blue", label = "Potential Energy")
    plt_E = plot(sol.t, round.(E_tot, sigdigits = 5), color = "darkgreen", title = "Total Energy", legend = false, titlefontsize = 10)

    # ===== Check if total energy has overcome the barrier
    lb, ub = extrema(E_tot)
    lb < 0 < ub ? push!(Aw_set, [A, ω]) : nothing
    # if ub > 0
    #   break
    # end
    # ===== Potential Energy profile
    xs_max = 7
    xs = -xs_max:0.1:xs_max
    V(x) = 1 / 2 * c * x^2 + 1 / 4 * β * x^4
    plt_V = plot(xs, V.(xs), ylims = [-E_barrier, E_barrier], title = "Potential Surface (Barrier: $E_barrier)", legend = false, titlefontsize = 10)

    # ===== Final plot
    l = @layout [a{0.01h}; b c; d e; f]
    plt = plot(suptitle, plt_ds, plt_V, plt_Ek, plt_Ep, plt_E,
      size = (800, 800),
      layout = l,
      # title = "Dufﬁng Oscillator",
      dpi = 600)
    # display(plt)
  end
end

plot([x[2] for x in record_fpt],[ x[1] for x in record_fpt], lw =3)
## Plot A w curve
MM = hcat(Aw_set...)'
MA = MM[:, 1]
Mw = MM[:, 2]
scatter(MA, Mw,
  xlabel = "switching Amplitude", ylabel = "switching frequency")





##
df = DataFrame(A = MA, B = Mw)

using CSV, DataFrames
# load data
df = CSV.File("./SciML_Notch/Duffing_Aw_c = -2, β = 1, A = 0.01:0.01:20, ω = 0.001:0.1:10_γ = 0.1.csv") |> DataFrame

using Plots
scatter(df.A, df.B,
  xlabel = "switching Amplitude", ylabel = "switching frequency")



















######### =================================================================
######### ====================== asymmetric potenial =======================
######### =================================================================
######### =================================================================
##
##
γ = 0.0001
c = -10
β = 1
A = 8.004
ω = 0.8
# θ = 3.0
θ = 4.0

# In an asymmetric potenial
function duffing_asy(du, u, p, t)
  x, v = u
  c, β, γ, A, ω, θ = p
  du[1] = u[2]
  du[2] = -c * u[1] - β * u[1]^3 - γ * u[2] - θ + A * cos(ω * t)
end
# V(x) = 1 / 2 * c * x^2 + 1 / 4 * β * x^4;    # potential function (symmetric)
V(x) = 1 / 2 * c * x^2 + θ * x + 1 / 4 * β * x^4  # potential function (asymmetric)
Vp(x) = c * x + θ + β * x^3
rts = values(IntervalRootFinding.roots(x -> Vp(x), 0 .. 10))
midpoints = mid.(interval.(rts))
# x_pos1 = 0.302775
# x_pos2 = 2.99999 # result from above line
@show x_pos2, x_pos1 = midpoints
E = V(x_pos1) - V(x_pos2)

u0 = [x_pos2; 0.0]
p = [c, β, γ, A, ω, θ]
Δt = 1e-2
T = 1000
tspan = (0.0, T)
prob = ODEProblem(duffing_asy, u0, tspan, p)
sol = solve(prob, RK4(), adaptive = false, dt = Δt)


# ====== Visualization
plt_ds = plot(sol, label = ["Position" "Velocity"],
  xlabel = "Time (t)", color = ["blue" "orange"], linewidth = [2.5 1.5],
  title = "Dynamics (A: $A, ω = $ω)", titlefontsize = 10)

pos = sol[1, :]
velocity = sol[2, :]
E_k = 1 / 2 * velocity .^ 2
E_p = V.(pos)
E_tot = 1 / 2 * velocity .^ 2 .+ V.(pos)

plt_Ek = plot(sol.t, E_k, color = "orange", label = "Kinetic Energy : " * L"E_{k}", xlabel = "Time (t)", ylabel = L"E_{k}")
plt_Ep = plot(sol.t, E_p, color = "blue", label = "Potential Energy" * L"E_{p}", xlabel = "Time (t)", ylabel = L"E_{p}")
plt_E = plot(sol.t, round.(E_tot, sigdigits = 5), color = "darkgreen", xlabel = "Time (t)", ylabel = L"E_{total}",
  label = "Total Energy", legend = true, titlefontsize = 10)

# ======  potential surface
xs = -5:0.1:5
plt_V = plot(xs, V.(xs), label = "Potential Surface",
  xlabel = "Position", dpi = 300)

l = @layout [b c; d e; f]
plt = plot(plt_ds, plt_V, plt_Ek, plt_Ep, plt_E,
  size = (800, 800),
  layout = l,
  # title = "Dufﬁng Oscillator",
  dpi = 600)
display(plt)
# savefig(plt, "./figures/asyn_duffing_sim.png")


stable_start = trunc(Int64, length(pos) * 3 / 4)
stable_pos = pos[stable_start:end]
dn, up = extrema(stable_pos)
(dn + up) / 2 < x_pos1

## ==============================================================================
## ==============================================================================
E_middle = V(x_pos1)
Aw_set = []
# suptitle = plot(title = "Duffing Oscillator", framestyle = nothing, showaxis = false, xticks = false, yticks = false, margin = 0mm)
# c = -10; β = 1; γ = 0.1; θ = 10.0;
@showprogress for ω = 0:0.5:15#collect(exp10.(-2:0.02:1))#0.01:0.5:30
  for A = 0:2:100#collect(exp10.(-2:0.02:1.4))
    p[1] = c
    p[2] = β
    p[3] = γ
    p[4] = A
    p[5] = ω
    p[6] = θ
    # @show A, ω
    prob1 = remake(prob, p = p)
    sol = solve(prob1, RK4(), adaptive = false, dt = Δt)

    ### ==== Visualization ======

    # ===== Dynamics plot
    # plt_ds = plot(sol, label = ["Position" "Velocity"],
    #   xlabel = "Time (t)",
    #   title = "Dynamics (A: $A, ω = $ω)", titlefontsize = 10)

    # # ===== Total energy plot
    pos = sol[1, :]
    velocity = sol[2, :]
    E_tot = 1 / 2 * velocity .^ 2 .+ V.(pos)

    # # check for steady states E
    # nc = 3 * Int(trunc(size(E_tot)[1] / 4))
    # E_tot_trunc = E_tot[nc:end]
    # plt_E = plot(sol.t, round.(E_tot, sigdigits = 5), title = "Total Energy", legend = false, titlefontsize = 10)

    # ===== Check if total energy has overcome the barrier
    # lb, ub = extrema(E_tot)
    # lb < E_middle < ub ? push!(Aw_set, [A, ω]) : nothing

    # lb, ub = extrema(E_tot_trunc)
    # E_middle < ub ? push!(Aw_set, [A, ω]) : nothing

    # ===== Criterion 1 for switching
    # pos_min = minimum(pos)
    # pos_min < x_pos1 ? push!(Aw_set, [A, ω]) : nothing
    #
    # if pos_min < x_pos1
    #     break
    # end

    # ===== Criterion 2 for switching
    stable_start = trunc(Int64, length(pos) * 3 / 4)
    stable_pos = pos[stable_start:end]
    dn, up = extrema(stable_pos)

    (dn + up) / 2 < x_pos1 ? push!(Aw_set, [A, ω]) : nothing

    if up < x_pos1
      break
    end
    # ===== Potential Energy profile
    # xs_max = 7
    # xs = -xs_max:0.1:xs_max
    # plt_V = plot(xs, V.(xs), title = "Potential Surface (Barrier: $E)", legend = false, titlefontsize = 10)

    # ===== Final plot
    # l = @layout [a{0.01h}; b; c d]
    # plt = plot(suptitle, plt_ds, plt_E, plt_V,
    #   layout = l,
    #   # title = "Dufﬁng Oscillator",
    #   dpi = 300)
    # display(plt)
  end
end




## Plot A w curve large region
MM = hcat(Aw_set...)'
MA = MM[:, 1]
Mw = MM[:, 2]
plt = scatter(MA, Mw,
  xlabel = "Switching Amplitude (A)", ylabel = "Switching frequency (ω)",
  m = (:viridis, 0.2, Plots.stroke(1, :green)),
  label = "swithing event", dpi = 500)
add_Boundary_event(MA, Mw, plt, lw = 1.5)
save_path = joinpath(pwd(), "figures", "A_w_collection/γ = 0.1, c = -10, β = 1, θ = 4.0/") # generate path to save
isdir(save_path) || mkpath(save_path)
savefig(plt, save_path * "Asym_A_w_region3.png")

## zoom region for low frequency
plt = scatter(MA, Mw,
  xlabel = "Switching Amplitude (A)", ylabel = "Switching frequency (ω)",
  xlims = [0, 100], ylims = [0, 15],
  m = (:viridis, 0.2, Plots.stroke(2, :green)),
  label = "swithing event", dpi = 500)
add_Boundary_event(MA, Mw, plt, lw = 1.5)
save_path = joinpath(pwd(), "figures", "A_w_collection/") # generate path to save
isdir(save_path) || mkpath(save_path)
savefig(plt, save_path * "Asym_A_w_zoomed_region.png")








## ================== Phase indenpent A-ϕ curve =================
γ = 1.1
c = -10
β = 1
A = 8.004
ω = 0.8
θ = 4.0
ϕ = 0.0

# In an asymmetric potenial
function duffing_asy(du, u, p, t)
  x, v = u
  c, β, γ, A, ω, θ, ϕ = p
  du[1] = u[2]
  du[2] = -c * u[1] - β * u[1]^3 - γ * u[2] - θ + A * cos(ω * t + ϕ)
end

V(x) = 1 / 2 * c * x^2 + θ * x + 1 / 4 * β * x^4  # potential function (asymmetric)
Vp(x) = c * x + θ + β * x^3
rts = values(IntervalRootFinding.roots(x -> Vp(x), 0 .. 10))
midpoints = mid.(interval.(rts))
@show x_pos2, x_pos1 = midpoints
E = V(x_pos1) - V(x_pos2)

u0 = [x_pos2; 0.0]
p = [c, β, γ, A, ω, θ, ϕ]
Δt = 1e-2
T = 1000
tspan = (0.0, T)
prob = ODEProblem(duffing_asy, u0, tspan, p)
sol = solve(prob, RK4(), adaptive = false, dt = Δt)
plot(sol)

##
Aw_set = []
function switched_late_stage_mean(pos)
  stable_start = trunc(Int64, length(pos) * 3 / 4)
  stable_pos = pos[stable_start:end]
  dn, up = extrema(stable_pos)
  average_osci = (dn + up) / 2
end
@showprogress for ω = 0:0.5:15#collect(exp10.(-2:0.02:1))#0.01:0.5:30
  for A = 0:2:100#collect(exp10.(-2:0.02:1.4))
    switch_set = []
    for ϕ in 0:2π
      p[1] = c
      p[2] = β
      p[3] = γ
      p[4] = A
      p[5] = ω
      p[6] = θ
      p[7] = ϕ

      prob1 = remake(prob, p = p)
      sol = solve(prob1, RK4(), adaptive = false, dt = Δt)

      ### ==== Visualization ======

      # ===== Dynamics plot
      # plt_ds = plot(sol, label = ["Position" "Velocity"],
      #   xlabel = "Time (t)",
      #   title = "Dynamics (A: $A, ω = $ω)", titlefontsize = 10)

      # # ===== Total energy plot
      pos = sol[1, :]
      velocity = sol[2, :]
      E_tot = 1 / 2 * velocity .^ 2 .+ V.(pos)

      # # check for steady states E
      # nc = 3 * Int(trunc(size(E_tot)[1] / 4))
      # E_tot_trunc = E_tot[nc:end]
      # plt_E = plot(sol.t, round.(E_tot, sigdigits = 5), title = "Total Energy", legend = false, titlefontsize = 10)

      # ===== Check if total energy has overcome the barrier
      # lb, ub = extrema(E_tot)
      # lb < E_middle < ub ? push!(Aw_set, [A, ω]) : nothing

      # lb, ub = extrema(E_tot_trunc)
      # E_middle < ub ? push!(Aw_set, [A, ω]) : nothing

      # ===== Criterion 1 for switching
      # pos_min = minimum(pos)
      # pos_min < x_pos1 ? push!(Aw_set, [A, ω]) : nothing
      #
      # if pos_min < x_pos1
      #     break
      # end

      # ===== Criterion 2 for switching
      # stable_start = trunc(Int64, length(pos)*3/4)
      # stable_pos = pos[stable_start:end]
      # dn, up = extrema(stable_pos)
      # average_osci = (dn+up)/2
      # average_osci < x_pos1 ? push!(Aw_set, [A, ω]) : nothing
      average_osci = switched_late_stage_mean(pos)
      average_osci < x_pos1 ? push!(switch_set, -1) : push!(switch_set, 1)
      if average_osci < x_pos1
        break
      end
      # ===== Potential Energy profile
      # xs_max = 7
      # xs = -xs_max:0.1:xs_max
      # plt_V = plot(xs, V.(xs), title = "Potential Surface (Barrier: $E)", legend = false, titlefontsize = 10)

      # ===== Final plot
      # l = @layout [a{0.01h}; b; c d]
      # plt = plot(suptitle, plt_ds, plt_E, plt_V,
      #   layout = l,
      #   # title = "Dufﬁng Oscillator",
      #   dpi = 300)
      # display(plt)
    end
    if -1 in switch_set
      push!(Aw_set, [A, ω])
    end

    if -1 in switch_set
      break
    end

  end
end
## Plot A w curve large region
MM = hcat(Aw_set...)'
MA = MM[:, 1]
Mw = MM[:, 2]
plt = scatter(MA, Mw,
  xlabel = "Switching Amplitude (A)", ylabel = "Switching frequency (ω)",
  m = (:viridis, 0.2, Plots.stroke(1, :green)), title = "γ : $γ",
  label = "swithing event", dpi = 500)
add_Boundary_event(MA, Mw, plt, lw = 1.5)
save_path = joinpath(pwd(), "figures", "A_w_collection/ϕ_independent/") # generate path to save
isdir(save_path) || mkpath(save_path)
savefig(plt, save_path * "Asym_A_w_γ : $γ.png")
