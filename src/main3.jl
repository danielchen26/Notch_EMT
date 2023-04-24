# This file is to explore the phase effects on epigenetic switching.
## ====================== Loading packages ==========================
using DifferentialEquations
using Plots
# using ModelingToolkit
using Catalyst
using ProgressMeter
using DataFrames, CSV, Random
using StatsPlots
using Latexify

## =====================  loading parameters data library =====================
db = CSV.File("../../../Notch_EMT_data/Notch_params_complete.csv") |> DataFrame

p_names = names(db)[1:11]
initi_names = names(db)[13:end]

df = db[shuffle(1:nrow(db))[1:10], :]
parameter_set = df[:, p_names]
initial_condition = df[:, initi_names]






## =====================  loading function library =====================
function make_cb(ts_in, index, value)
    ts = ts_in
    condition(u, t, integrator) = t in ts
    function affect!(integrator)
        if integrator.t == ts[1]
            integrator.p[index] = value
        elseif integrator.t == ts[2]
            integrator.p[index] = 0.0
        end
    end
    cb = DiscreteCallback(condition, affect!, save_positions = (true, true))
    return ts, cb
end

"""
test if the steady state switched
	-1 : switched
	 1 : not swithced
"""
function check_switching(sol)
    H4_i, H27_i = sol(ts[1])[[6, 9]]
    H4_f, H27_f = sol(tspan[2])[[6, 9]]
    init = H4_i / H27_i < 1 ? 1 : -1
    after = H4_f / H27_f < 1 ? 1 : -1
    return init * after
end

function add_Boundary_event(switch_amplitude, switch_frequency, plt)
    df = DataFrame(switch_amplitude = switch_amplitude, switch_frequency = switch_frequency)
    # CSV.write("switching_freq_amplitude.csv",df)

    gp = groupby(df, :switch_frequency)
    boundary_amplitude = []
    for i in gp
        append!(boundary_amplitude, minimum(i.switch_amplitude))
    end

    critical_freq = [keys(gp)[i].switch_frequency for i = 1:gp.ngroups]
    boundary_amplitude
    plt_boundary = plot!(plt, boundary_amplitude, critical_freq, label = "switching boundary", lw = 3)
end

# save A-w data
function save_A_w_data(switch_amplitude, switch_frequency, path)
    df = DataFrame(switch_amplitude = switch_amplitude, switch_frequency = switch_frequency)
    CSV.write(path * "switch_amplitude.csv", df)
end


# save initial condition, parameter set data
function save_init_param(parameter_set, initial_condition, id, path)
    CSV.write(path * "initial_condition.csv", DataFrame(initial_condition[id, :]))
    CSV.write(path * "parameter_set.csv", DataFrame(parameter_set[id, :]))
end






## ===================== Main Model ========================
model = @reaction_network begin
    # (A*(1+sign(cos(w*t))),1.0),     R ↔ NR               					# NICD binds RBPJ
    (A * (1 + cos(w * t + ϕ)), 1.0), R ↔ NR               # NICD binds RBPJ
    (k1, k2), M + R ↔ MR          # MITF binds RBPJ
    k0, MR --> MR + KDM5A            # MITF-RBPJ recruit KDM5A
    d, H4 + KDM5A --> H0 + KDM5A       # Demethylation of active mark
    m, H0 + PRC2 --> H27 + PRC2   # Methylation to get repressive mark
    1.0, H27 + KDM6A --> H0 + KDM6A  # Demethylation of repressive mark
    1.0, H0 + KMT --> H4 + KMT    # Methylation to get active mark
    p, H27 --> H27 + PRC2           # PRC2 is enhenced by H27 mark
    kk, H4 --> H4 + KDM6A        # KDM6A is enhenced by H4 mark
    pp, H4 --> H4 + KMT          # KMT is enhenced by H4 mark
    k, H27 --> H27 + KDM5A                # KDM5A is enhenced by H27 mark
    δ, (PRC2, KDM5A, KDM6A, KMT) --> ∅                    # Degradation of histone reader and writers
    α1, ∅ --> (KDM6A, KMT, PRC2, KDM5A)
end k0 k1 k2 d m p k pp kk δ α1 w A ϕ # put A at last as the control 13th variable
@show species(model)
@show parameters(model)
@show speciesmap(model)
@show paramsmap(model)
@show latexify(model)
# ODE_equations = convert(ODESystem, model)
# @show latexify(ODE_equations)


## ================================ Single Run for the model ============================================
# example for single run with phase
# 11-element Vector{Pair{Term{Real, Base.ImmutableDict{DataType, Any}}, Float64}}:
#      R(t) => 6.588723439378914
#     NR(t) => 8.472870571384637e-18
#      M(t) => 6.588723439378912
#     MR(t) => 43.41127656062109
#  KDM5A(t) => 44.41127656062109
#     H4(t) => 0.006677187301257658
#     H0(t) => 0.11326330551407315
#   PRC2(t) => 799.0809521149547
#    H27(t) => 49.88005950718467
#  KDM6A(t) => 1.1068349968201225
#    KMT(t) => 1.0734490603138342

# 14-element Vector{Pair{Sym{Real, Base.ImmutableDict{DataType, Any}}, Float64}}:
#  k0 => 1.0
#  k1 => 1.0
#  k2 => 1.0
#   d => 0.41
#   m => 0.61
#   p => 16.0
#   k => 0.0
#  pp => 11.0
#  kk => 16.0
#   δ => 1.0
#  α1 => 1.0
#   w => 0.3
#   A => 0.0
#   ϕ => 2.0943951023931953
u0 = collect(initial_condition[1, :])
phase = 0#2π / 3
p = vcat([collect(parameter_set[1, :]), 0.1, 0, phase]...)
ts, cb = make_cb([50.0, 70.0], 13, 800060)
tspan = (0.0, 350.0)
u0map = species(model) .=> u0
pmap = parameters(model) .=> p
prob = ODEProblem(model, u0map, tspan, pmap)
@time sol = solve(prob, Rosenbrock23(), callback = cb, tstops = ts)
plot(sol,
    vars = [4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
    lw = 1.5,
    xlabel = "Time", ylabel = "Concentration",
    title = "Switching Dynamics",
    dpi = 500)



## ==== amination for phase effects =================================
# note: this function is following the above example.
function phase_animation()
    anim = @animate for phase ∈ 0:0.2:2π
        phase = round(phase, digits = 2)
        # u0 = collect(initial_condition[1, :])
        p = vcat([collect(parameter_set[1, :]), 0.3, 0, phase]...)
        # ts, cb = make_cb([50.0, 100.0], 13, 300.0)
        # tspan = (0.0, 150.0)
        # u0map = species(model) .=> u0
        pmap = parameters(model) .=> p
        prob = ODEProblem(model, u0map, tspan, pmap)
        sol = solve(prob, Rosenbrock23(), callback = cb, tstops = ts)
        plot(sol,
            vars = [4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
            lw = 1.5,
            xlabel = "Time", ylabel = "Concentration",
            title = "Switching Dynamics with phase : $phase",
            dpi = 500)
    end
    save_path = joinpath(@__DIR__, "figures", "phase_animation/") # generate path to save
    isdir(save_path) || mkpath(save_path) # create directory if not exist
    CSV.write(save_path * "initial_condition.csv", DataFrame(var = initi_names, values = u0))
    CSV.write(save_path * "parameter_set.csv", DataFrame(var = vcat([p_names, "ω", "A", "ϕ"]...), values = p))
    gif(anim, save_path * "phase_anim_critical.gif", fps = 3)
end

phase_animation()





##  ================================ With random parameters generate switching dynamics ============================================
# df = db[shuffle(1:nrow(db))[1:10], :]
parameter_set = db[:, p_names]
initial_condition = db[:, initi_names]
rand_idx_set = shuffle(1:nrow(db))
@show db_idx = rand(rand_idx_set)
@show parameter_set[db_idx, :]

## ===========================================================================================
# example for single run with phase
freq = 0.1
amplitude = 50.0
p = vcat([collect(parameter_set[db_idx, :]), freq, 0.0, 0.0]...)
u0 = collect(db[db_idx, initi_names])
ts, cb = make_cb([50.0, 100.0], 13, amplitude)
tspan = (0.0, 150.0)
u0map = species(model) .=> u0
pmap = parameters(model) .=> p
prob = ODEProblem(model, u0map, tspan, pmap)
@time sol = solve(prob, Rosenbrock23(), callback = cb, tstops = ts)
plt = plot(sol,
    # vars = [4, 6, 9, 5, 10], # [MR,KDM5A,H4,H27,KDM6A]
    vars = [4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
    lw = 2.5,
    xlabel = "Time", ylabel = "Concentration",
    foreground_color_legend = nothing,
    # title = "Switching Dynamics",
    dpi = 300)
plot!(plt, [0, ts[1], ts[2], tspan[end]], [0, 0, amplitude, 0],
    label = "Input", seriestype = :steppre, line = (:dashdot, 2), alpha = 0.8,
    fill = (0, 0.2, :purple), color = "black", dpi = 300)

## ===== save plot and generate paths
class = "activation_time/"
# class = "oscillating_input1/" # saved two folder contains 2 ω near the transition.
save_path = joinpath(@__DIR__, "figures", "switching_dynamics/$class") # generate path to save
isdir(save_path) || mkpath(save_path) # create directory if not exist
CSV.write(save_path * "initial_condition.csv", DataFrame(var = initi_names, values = u0))
CSV.write(save_path * "parameter_set.csv", DataFrame(var = vcat([p_names, "ω", "A", "ϕ"]...), values = p))
savefig(plt, save_path * "switching_dynamics_const_A=50_w=0.1.png")





## =============== A-w-ϕ curve ================
switch_amplitude = []
switch_frequency = []
switch_phase = []
tspan = [0.0, 150.0]
@showprogress for phase = 0:2π
    for freq in exp10.(-4:0.1:0)
        for amplitude = 0:10:1000#10:0.05:30
            # frequency modulation
            # p1 = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, freq, 0.0]

            # p = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, freq, 0.0, phase]
            # u0 = [6.0, 0.0, 6.0, 40.0, 500.0, 0.6876, 0.678, 500.0, 50.6344, 1.0, 2.0]

            # u0 = [6.59, 0.0, 6.59, 43.41, 61.47, 0.02, 0.54, 48.17, 49.43, 1.09, 1.09]
            # @show parameter_set[db_idx, :]
            # @show initial_condition[db_idx, :]

            p = vcat([collect(parameter_set[db_idx, :]), freq, 0.0, phase]...)
            pmap = parameters(model) .=> p
            u0 = collect(initial_condition[db_idx, :])
            u0map = species(model) .=> u0
            ts, cb = make_cb([50.0, 100.0], 13, amplitude)
            prob1 = ODEProblem(model, u0map, tspan, pmap)
            sol1 = solve(prob1, Rosenbrock23(), callback = cb, tstops = ts)

            # plt1 = plot(sol1, vars = [5, 6, 9, 10], lw = 1.5, title = "Amp: $amplitude, Frequency: $freq")
            # display(plt1)

            check = check_switching(sol1)
            if check == -1
                append!(switch_amplitude, amplitude)
                append!(switch_frequency, freq)
                append!(switch_phase, phase)
            end

            if check_switching(sol1) == -1
                break
            end
        end
    end
end


plt = plot(switch_amplitude, switch_frequency, seriestype = :scatter, #yaxis=:log10,
    label = "switching events", title = "Frequency vs Amplitude at the switching",
    xlabel = "switch_amplitude", ylabel = "switch_frequency", dpi = 500)
add_Boundary_event(switch_amplitude, switch_frequency, plt)



##
plot(switch_amplitude, switch_frequency, switch_phase, zcolor = reverse(switch_phase), m = (4, 2, :starrynight, Plots.stroke(2)),
    xlabel = "A", ylabel = "ω", zlabel = "ϕ", title = "A-w-ϕ curve",
    dpi = 500)


df3d = DataFrame(freq = switch_frequency, phase = switch_phase, amp = switch_amplitude)
using VegaLite
df3d |>
@vlplot(
    :point,
    x = :amp,
    y = {"freq:q", order = "descending"},
    color = {"phase:o", scale = {scheme = "blueorange"}},
    width = 400,
    height = 400,
    # row = :phase
)


df3d_sub1 = filter(row -> row.phase in [0], df3d)
df3d_sub2 = filter(row -> row.phase in [1, 2, 3], df3d)
df3d_sub21 = filter(row -> row.phase in [0, 1, 2, 3], df3d)
df3d_sub3 = filter(row -> row.phase in [4, 5, 6], df3d)


#  ==== Plots A-w curves for all phase ====
plt = @df df3d plot(
    :amp,
    :freq,
    group = :phase,
    palette = :RdYlGn_7,
    m = (2, 4),
    legend_title = "phase",
    legend_position = :outertopright,
    # xlabel = "Switching Amplitude",
    # ylabel = "Switching Frequency"
    dpi = 500,
    # bg = RGB(0.2, 0.2, 0.5)
)
xlabel!(plt, "Switching Amplitude")
ylabel!(plt, "Switching Frequency")
savefig(plt, "phase_shift.png")

## === plot A-w curves for separate phase ====
#
# scheme = cgrad(:RdYlGn_11, 11, categorical = true)
plt1 = @df df3d_sub1 plot(
    :amp,
    :freq,
    group = :phase,
    # palette = scheme,
    c = ["#007cd4"],
    m = (2, 4),
    legend_title = "phase",
    legend_position = :topleft,
    # xlims = [0, 30],
    xlabel = "Switching Amplitude",
    ylabel = "Switching Frequency",
    foreground_color_legend = nothing,
    dpi = 500,
    # bg = RGB(0.2, 0.2, 0.5)
)
savefig(plt1, "./figures/A-w_phase0.png")
plt2 = @df df3d_sub2 plot(
    :amp,
    :freq,
    group = :phase,
    # palette = scheme,
    color = ["#ff6232" "#ffaa4f" "#ffdf7d"],
    m = (2, 4),
    legend_title = "phase",
    legend_position = :topleft,
    # xlims = [0, 30],
    xlabel = "Switching Amplitude",
    ylabel = "Switching Frequency",
    foreground_color_legend = nothing,
    dpi = 500,
    # bg = RGB(0.2, 0.2, 0.5)
)
savefig(plt2, "./figures/A-w_phase1-3.png")
#007cd4
plt21 = @df df3d_sub21 plot(
    :amp,
    :freq,
    group = :phase,
    # palette = scheme,
    color = ["#007cd4" "#6541b2" "#ba3a4c" "#b38400"],
    m = (2, 4),
    legend_title = "phase",
    legend_position = :topleft,
    # xlims = [0, 30],
    xlabel = "Switching Amplitude",
    ylabel = "Switching Frequency",
    foreground_color_legend = nothing,
    dpi = 500,
    # bg = RGB(0.2, 0.2, 0.5)
)
savefig(plt21, "./figures/A-w_phase0-3.png")
plt3 = @df df3d_sub3 plot(
    :amp,
    :freq,
    group = :phase,
    # palette = scheme,
    color = ["#97db57" "#3ec057" "#006a31"],
    m = (2, 4),
    legendtitle = "phase",
    legend_position = :topleft,
    # xlims = [0, 30],
    xlabel = "Switching Amplitude",
    ylabel = "Switching Frequency",
    foreground_color_legend = nothing,
    dpi = 500,
    # bg = RGB(0.2, 0.2, 0.5)
)
savefig(plt3, "./figures/A-w_phase4-6.png")
# plt_combo = plot(plt1, plt2, plt3, layout = (1, 3), size = (1800, 600))
# savefig(plt_combo, "phase_shift_combo.png")





## ===## =============== A-w-ϕ curve if ϕ is not controllable================
switch_amplitude = []
switch_frequency = []
switch_phase = []
tspan = [0.0, 200.0]

@showprogress for freq in 0:0.001:0.5#exp10.(-4:0.05:1)
    for amplitude = 0:0.2:50#15:0.1:20#14:0.01:18#10:0.05:30
        switch_set = []
        for phase = 0:2π
            # frequency modulation
            # p1 = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, freq, 0.0]
            # p = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, freq, 0.0, phase]
            # u0 = [6.0, 0.0, 6.0, 40.0, 500.0, 0.6876, 0.678, 500.0, 50.6344, 1.0, 2.0]
            # u0 = [6.59, 0.0, 6.59, 43.41, 61.47, 0.02, 0.54, 48.17, 49.43, 1.09, 1.09]
            # @show parameter_set[db_idx, :]
            # @show initial_condition[db_idx, :]

            p = vcat([collect(parameter_set[db_idx, :]), freq, 0.0, phase]...)
            pmap = parameters(model) .=> p
            u0 = collect(initial_condition[db_idx, :])
            u0map = species(model) .=> u0
            ts, cb = make_cb([50.0, 90.0], 13, amplitude)
            prob1 = ODEProblem(model, u0map, tspan, pmap)
            sol1 = solve(prob1, Rosenbrock23(), callback = cb, tstops = ts)

            # plt1 = plot(sol1, vars = [5, 6, 9, 10], lw = 1.5, title = "Amp: $amplitude, Frequency: $freq")
            # display(plt1)

            check = check_switching(sol1)
            append!(switch_set, check)
            # if check == -1 # if there exist a phase
            #     append!(switch_amplitude, amplitude)
            #     append!(switch_frequency, freq)
            #     # append!(switch_phase, phase)
            # end
            if check_switching(sol1) == -1
                break
            end
        end
        # @show switch_set
        if -1 in switch_set
            append!(switch_amplitude, amplitude)
            append!(switch_frequency, freq)
        end
        if -1 in switch_set
            break
        end
    end
end



##
plt = plot(switch_amplitude, switch_frequency, seriestype = :scatter, #yaxis=:log10,
    label = "switching events",# title = "Frequency vs Amplitude at the switching",
    xlabel = "Switch Amplitude", ylabel = "Switch Frequency", dpi = 500)
add_Boundary_event(switch_amplitude, switch_frequency, plt)
savefig(plt, "./figures/A-w_curve_ϕ_indenpendent_ID:1233_T=40.png")
