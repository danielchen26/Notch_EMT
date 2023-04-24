# This file is to explore the pulsatile effects on epigenetic switching.
## ====================== Loading packages and data library==========================
using Revise
using DifferentialEquations
using Plots;
gr(fontfamily="Helvetica");
using Catalyst
using Catalyst: parameters

using ProgressMeter
using DataFrames, CSV, Random, Distributions
using StatsPlots
using Latexify, Measures, FLoops, LaTeXStrings
includet("./Functions.jl")
db, p_names, initi_names, parameter_set, initial_condition, rand_idx_set, db_idx = loading_database()


## ===================== Main Model ========================

# model = @reaction_network begin 
#     # (A * (1 + sign(cos(w * t + ϕ))), 1.0), R ↔ NR               # NICD binds RBPJ
#     (A * (abs(cos(w * t + ϕ))), 1.0), R ↔ NR               # NICD binds RBPJ 🍏🔴
#     (k1, k2), M + R ↔ MR          # MITF binds RBPJ
#     k0, MR --> MR + KDM5A            # MITF-RBPJ recruit KDM5A
#     d, H4 + KDM5A --> H0 + KDM5A       # Demethylation of active mark
#     m, H0 + PRC2 --> H27 + PRC2   # Methylation to get repressive mark
#     1.0, H27 + KDM6A --> H0 + KDM6A  # Demethylation of repressive mark
#     1.0, H0 + KMT --> H4 + KMT    # Methylation to get active mark
#     p, H27 --> H27 + PRC2           # PRC2 is enhenced by H27 mark
#     kk, H4 --> H4 + KDM6A        # KDM6A is enhenced by H4 mark
#     pp, H4 --> H4 + KMT          # KMT is enhenced by H4 mark
#     k, H27 --> H27 + KDM5A                # KDM5A is enhenced by H27 mark
#     δ, (PRC2, KDM5A, KDM6A, KMT) --> ∅                    # Degradation of histone reader and writers
#     α1, ∅ --> (KDM6A, KMT, PRC2, KDM5A)
# end k0 k1 k2 d m p k pp kk δ α1 w A ϕ # put A at last as the control 13th variable
# @show species(model)
# @show parameters(model)
# @show speciesmap(model)
# @show paramsmap(model)

# @show latexify(model)
# ODE_equations = convert(ODESystem, model)
# @show latexify(ODE_equations)
# Graph(model)



## ============================= Import models =================================
signal_type = "bump"
model_bump = Import_model(; type=signal_type)
@show equations(model_bump)
signal_type = "pulsatile"
model_pulsatile = Import_model(; type=signal_type)
@show equations(model_pulsatile)

rn_latex, ode_latex = ode_latex_gen(model_pulsatile)

## ============================== single run for the model =============================
# 🔴
db_idx = 592
freq = 0.0;
phase = 0;
amplitude = 220;
T_init = 0.001;
ΔT = 100;
tspan = (0.0, 350.0);
@show db[db_idx, p_names]
model = model_pulsatile

u0map, pmap, p, tspan, ts, cb, sol = single_solve(; db_idx=db_idx, freq=freq, phase=phase, amplitude=amplitude, T_init=T_init, ΔT=ΔT, tspan=tspan, phase_reset=true)
plt = plot(sol, vars=[4, 5, 6, 9], lw=1.5, xlabel="Time", ylabel="Concentration", dpi=500)


## ====== I want to find if 🔴 changing the prc2 kinetic rate is able to result in early activation
prob_new = remake_prob(model, u0map, tspan, p; prc2=0.59)
@time sol = solve(prob_new, Rosenbrock23(), callback=cb, tstops=ts)
plt = plot(sol,
    vars=[4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
    lw=1.5,
    xlabel="Time", ylabel="Concentration",
    title="PRC2 rate : 0.4",
    dpi=500)
anim_prc2 = anim_prc2_changing(0:0.1:0.8, tspan=[0, 350], u0=u0map)



## ===========================================================================================
# --------- Single run with phase unreseted #? This is just a test
phase = 0
freq = 0.2
amplitude = 298
T_init = 100
ΔT = 10^2
tspan = (0.0, 3 * ΔT)
db_idx = 592 # gene 1 with Dll4 const A, no matter how long the signal is giving to the system, no switch
# db_idx =49 # gene 2 with Dll1

# test---- solve step by step
p = vcat([collect(parameter_set[db_idx, :]), freq, 0.0, phase]...)
model = model_pulsatile
pmap = parameters(model) .=> p
u0 = collect(initial_condition[db_idx, :])
u0map = species(model) .=> u0
ts, cb = make_cb([T_init, T_init + ΔT], 13, amplitude)
prob1 = ODEProblem(model, u0map, tspan, pmap)
sol = solve(prob1, Rosenbrock23(), callback=cb, tstops=ts)
plot(sol)
## =================================================================







## ! skip this part ===========================================================================================
u0map, pmap, p, tspan, ts, cb, sol, phase = single_solve(; model=model_pulsatile, db_idx=db_idx, freq=freq, phase=phase, amplitude=amplitude, T_init=T_init, ΔT=ΔT, tspan=tspan);

@show t_switching = switching_time(; sol=sol, pulse_period=T_init:0.1:T_init+ΔT, idx=[6, 9], return_plot=false)

plt = plot(sol,
    vars=[4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
    lw=2,
    xlabel="Time", ylabel="Concentration",
    foreground_color_legend=nothing,
    # title = "A : $amplitude, freq = $freq, switch time : $t_switching",
    dpi=500)

# plot!(plt, [0, ts[1], ts[2], tspan[end]], [0, 0, amplitude, 0],
# label = "Sustainable Input", seriestype = :steppre, line = (:dashdot, 2), alpha = 0.8,
# # ylims = [0, 400],
# fill = (0, 0.3, :blue), color = "black", dpi = 300)

# # plot oscillating signal region
osci_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)))
# osci_signal(t, A, w, ϕ) = A * (abs(cos(w * t + ϕ))) #🔴
tt = ts[1]:0.01:ts[2]
osci_signal.(tt, amplitude, freq, 0.0)
plot!(plt, tt, osci_signal.(tt, amplitude, freq, 0.0),
    label="Pulsatile Input", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
    # ylims = [0, 700],
    fill=(0, 0.3, :darkgreen), color="black", dpi=300)

# save_plot(plt;filename = "switching_dynamics")


## ====== to wrap up the above code =================================
# function test(freq)
#         # --------- Single run with phase for paper plots
#     phase = 0
#     freq = freq
#     amplitude = 298
#     T_init = 100
#     ΔT = 10^2
#     tspan = (0.0, 3 * ΔT)
#     db_idx = 592 # gene 1 with Dll4 const A, no matter how long the signal is giving to the system, no switch
#     # db_idx =49 # gene 2 with Dll1

#     u0map, pmap, p, tspan, ts, cb, sol = single_solve(; model=model_pulsatile, db_idx=db_idx, freq=freq, phase=phase, amplitude=amplitude, T_init=T_init, ΔT=ΔT, tspan=tspan)

#     @show t_switching = switching_time(; sol=sol, pulse_period=T_init:0.1:T_init+ΔT, idx=[6, 9], return_plot=false)

#     plt = plot(sol,
#         vars=[4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
#         lw=2,
#         xlabel="Time", ylabel="Concentration",
#         foreground_color_legend=nothing,
#         # title = "A : $amplitude, freq = $freq, switch time : $t_switching",
#         dpi=500)

#     # plot!(plt, [0, ts[1], ts[2], tspan[end]], [0, 0, amplitude, 0],
#     # label = "Sustainable Input", seriestype = :steppre, line = (:dashdot, 2), alpha = 0.8,
#     # # ylims = [0, 400],
#     # fill = (0, 0.3, :blue), color = "black", dpi = 300)

#     # # plot oscillating signal region
#     osci_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)))
#     # osci_signal(t, A, w, ϕ) = A * (abs(cos(w * t + ϕ))) #🔴
#     tt = ts[1]:0.01:ts[2]
#     osci_signal.(tt, amplitude, freq, 0.0)
#     plot!(plt, tt, osci_signal.(tt, amplitude, freq, 0.0),
#         label="Pulsatile Input", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
#         # ylims = [0, 700],
#         fill=(0, 0.3, :darkgreen), color="black", dpi=300)
#     display(plt)
# end

# for freq in 0:0.01:0.5
#     test(freq)
# end




## ====== Single plot for 1 gene comparing pulsatile case vs bump case=======
T_init = 1e-10
# ======= pulsatile model
model = model_pulsatile
signal_type = "pulsatile"
plt_pulsatile = single_solve_plot(; model=model, db_idx=49, phase=0, freq=0.5, amplitude=65.0, T_init=T_init, ΔT=80, type=signal_type)
# ======== bump model
model = model_bump
signal_type = "bump"
plt_bump = single_solve_plot(; model=model, db_idx=49, phase=0, freq=0.5, amplitude=65.0, T_init=T_init, ΔT=80, type=signal_type)
plt_2model = plot(plt_pulsatile, plt_bump, layout=(2, 1))
##


# write a function to plot 1 gene comparing pulsatile case vs bump case
function single_solve_plot_pulsatile_bump(; db_idx=49, phase=0, freq=0.5, amplitude=65.0, T_init=0.01, ΔT=100)
    plt_pulsatile = single_solve_plot(; model=model_pulsatile, db_idx=db_idx, phase=phase, freq=freq, amplitude=amplitude, T_init=T_init, ΔT=ΔT, type="pulsatile")
    plt_bump = single_solve_plot(; model=model_bump, db_idx=db_idx, phase=phase, freq=freq, amplitude=amplitude, T_init=T_init, ΔT=ΔT, type="bump")
    plot(plt_pulsatile, plt_bump, layout=(2, 1))
end


#! Have to set T_init to 0.01 to avoid the discontinuity
for ΔT ∈ 0.01:10:1000
    plt = single_solve_plot_pulsatile_bump(ΔT=ΔT, T_init=0.01)
    display(plt)
end





## ====== single plot for 1 gene ID 592 ======= Dll4 vs Dll1 within the first the duartiona of pulse given T_init
T_init = 1e-10
# ======= sustained model
plt_sustained = single_solve_plot(; model=model_pulsatile, db_idx=592, phase=0, freq=0.0, amplitude=165.0, T_init=T_init, ΔT=50, type="sustained", phase_reset=true)
# ======= pulsatile model
plt_pulsatile = single_solve_plot(; model=model_pulsatile, db_idx=592, phase=0, freq=0.010, amplitude=165.0, T_init=T_init, ΔT=50, type="pulsatile", phase_reset=true)

plot(plt_sustained, plt_pulsatile, layout=(2, 1))
##


# * write function for ploting 1 gene ID default to 592 comparing sustained signal vs pulsatile signal.
function single_solve_plot_sustained_pulsatile(; db_idx=592, phase=0, freq=0.0, amplitude=165.0, T_init=0.01, ΔT=100)
    plt_sustained = single_solve_plot(; model=model_pulsatile, db_idx=db_idx, phase=0.0, freq=0.0, amplitude=amplitude, T_init=T_init, ΔT=ΔT, type="sustained", phase_reset=true)
    plt_pulsatile = single_solve_plot(; model=model_pulsatile, db_idx=db_idx, phase=phase, freq=freq, amplitude=amplitude, T_init=T_init, ΔT=ΔT, type="pulsatile", phase_reset=true)
    plot(plt_sustained, plt_pulsatile, layout=(2, 1))
end


# * showed some freq and amplitude allows pulsatile signal switch states after signal is off (between pulses)
for freq ∈ 0.01:0.05:1#, amplitude ∈ 165:10:300
    plt = single_solve_plot_sustained_pulsatile(ΔT=100, T_init=0.01, freq=freq)#, amplitude = amplitude)
    display(plt)
    sleep(0.1)
end





## ======================== to find index with A = 64, Dll4 switch 
for id = 100:200
    find_id_Dll4_vs_Dll1(id, amplitude=300, prc2=0.64)
    sleep(0.1)
end
## ======================== to find index with A = 64, Dll4 switch





## * Exmaple 1 ======= Dll4 vs Dll1 for gene id:49 ====
single_gene_id = 49
id2_freq = 0.15
# phase2 = 5
amplitude1 = 62
amplitude2 = 62
T_init = 1e-10
ΔT = 100
prc2 = 0.64
plotly()
plt_gene1_Dll4, plt_gene1_Dll1, plt_2genes_compare_id_49 =
    Two_Genes_TS_by_Prc2(;
        # model = model_bump,
        model=model_pulsatile,
        id1=single_gene_id, id2=single_gene_id,
        id2_freq=id2_freq, phase2=phase2, amplitude1=amplitude1,
        amplitude2=amplitude2, prc2=prc2,
        T_init=T_init, ΔT=ΔT, title_on=true, legend_title_on=false,
        vars_to_show=[5, 6, 9], #tspan_rt = 2, # 
        type="pulsatile",
        # type="bump",
        phase2_reset=true) #🔴 specify the type
plt_gene1_Dll4
plt_gene1_Dll1
plt_2genes_compare_id_49
# savefig(plt_2genes_compare_id_49,"./figures/APS_id_49_Dll4_Dll1_compare_bump.png")
##
# savepath = pathgen("pulsatile")
# savefig(plt_gene1_Dll4, savepath * "plt2_gene1_Dll4.png")
# savefig(plt_gene1_Dll1, savepath * "plt2_gene1_Dll1.png")




## * Example 2 ======= Dll4 vs Dll1 for gene id:592 ====
single_gene_id = 592
id2_freq = 0.09
amplitude1 = amplitude2 = 134
T_init = 1e-10
ΔT = 100
prc2 = 0.41

plt_gene2_Dll4, plt_gene2_Dll1, plt_2genes_compare_id_592 =
    Two_Genes_TS_by_Prc2(; model=model_pulsatile,
        id1=single_gene_id, id2=single_gene_id,
        id2_freq=id2_freq, phase2=phase2, amplitude1=amplitude1,
        amplitude2=amplitude2, prc2=prc2,
        T_init=T_init, ΔT=ΔT, title_on=true, legend_title_on=true,
        type="pulsatile",
        phase2_reset=true) #🔴 specify the type
plt_gene2_Dll4
plt_gene2_Dll1
plt_2genes_compare_id_592
##
# savepath = pathgen("bump")
# savefig(plt_gene2_Dll4, savepath * "plt_gene2_Dll4.png")
# savefig(plt_gene2_Dll1, savepath * "plt_gene2_Dll1.png")






## *======= Animation for gene ID:49 Dll4 vs Dll1, when varying PRC2 rates ===========  
gene_id = 49
id2_freq = 0.13
amplitude = 100
anim_comp = Two_genes_prc2_TS_animation(; prc2_range=0:0.02:1, model=model_pulsatile,
    id1=gene_id,
    id2=gene_id,
    id2_freq=id2_freq,
    amplitude1=amplitude,
    amplitude2=amplitude)
gif(anim_comp, fps=1)
##

## * gene ID: 49 as an example, visualization for the two genes TS by changing prc2 rate 

_, _, plt_2genes_compare = Two_Genes_TS_by_Prc2(; id1=49, id2=49, id2_freq=0.13, prc2=0.8, amplitude1=100, amplitude2=100)
plt_2genes_compare


##

# savefig(plt_2genes_compare, "./figures/592_49_prc2=0.862.png")


## ------- loop over each parameter set to find large TS that is close to the ΔT[end]
# phase = 0
# freq = 0.0
# amplitude = 120
# T_init = 100
# ΔT = 100
# tspan = (0.0, 10^6)
# db_idx = 592
# for ΔT in 1:10^4:10^6
#     @show db_idx
#     @show ΔT
#     u0map, pmap, p, tspan, ts, cb, sol = single_solve(; db_idx = db_idx, freq = freq, phase = phase, amplitude = amplitude, T_init = T_init, ΔT = ΔT, tspan = tspan)
#     t_switching = switching_time(; sol = sol, pulse_period = T_init:0.1:T_init+ΔT, idx = [6, 9], return_plot = false)
#     @show t_switching
#     plt = plot(sol,
#         vars = [4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
#         lw = 2,
#         xlabel = "Time", ylabel = "Concentration",
#         foreground_color_legend = nothing,
#         title = "A : $amplitude, freq = $freq, switch time : $t_switching, ID : $db_idx",
#         dpi = 300)
#     display(plt)
# end

## -----🔶 animation of fixed amplitude and phase with increase frequency, and its switching time ----------------
# ------🔹 phase dependent case with phase ϕ = 0
# anim_freq_tswitch(; range=0:0.02:0.4, amplitude=400, db_idx=db_idx)




## ================  Calculate the A ω ϕ st relation for a single gene. specify the db_idx number for a gene.========================
db_idx = 49
# db_idx = 600 # 🍏 paper figure 6
df = A_ω_st_relation(; model=model_pulsatile,
                    db_idx=db_idx,
                    amplitude_range=0:50:300, 
                    freq_range=0:0.02:2,
                    ΔT=100)
plt_freq_vs_ST = df_freq_vs_ST_groupby_amp(df; amplitude_select = [], palette=cgrad([:goldenrod1, :dodgerblue4, :chartreuse3]), save=false)



# ! not does need this, no phase is needed now
## =============== A-w-ϕ curve ================
switch_amplitude = []
switch_frequency = []
switch_phase = []
tspan = [0.0, 150.0]

@time @showprogress for phase = 0:2π
    for freq in exp10.(-4:0.1:0)
        for amplitude = 0:10:1500#10:0.05:30
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
            ts, cb = make_cb([50.0, 50.0 + ΔT], 13, amplitude)
            prob1 = ODEProblem(model, u0map, tspan, pmap)
            sol1 = solve(prob1, Rosenbrock23(), callback=cb, tstops=ts)

            # plt1 = plot(sol1, vars = [5, 6, 9, 10], lw = 1.5, title = "Amp: $amplitude, Frequency: $freq")
            # display(plt1)

            check = check_switching(sol1, ts, tspan)
            if check == -1
                append!(switch_amplitude, amplitude)
                append!(switch_frequency, freq)
                append!(switch_phase, phase)
            end

            if check_switching(sol1, ts, tspan) == -1
                break
            end
        end
    end
end


plt = plot(switch_amplitude, switch_frequency, seriestype=:scatter, #yaxis=:log10,
    label="switching events", title="Frequency vs Amplitude at the switching",
    xlabel="switch_amplitude", ylabel="switch_frequency", dpi=500)
add_Boundary_event(switch_amplitude, switch_frequency, plt)


plt = plot(df4d.amp, df4d.freq, seriestype=:scatter, #yaxis=:log10,
    label="switching events", title="Frequency vs Amplitude at the switching",
    xlabel="switch_amplitude", ylabel="switch_frequency", dpi=500)

add_Boundary_event(df4d.amp, df4d.freq, plt)

##
df3d = DataFrame(freq=switch_frequency, phase=switch_phase, amp=switch_amplitude)
df3d_sub1 = filter(row -> row.phase in [0], df3d)
df3d_sub2 = filter(row -> row.phase in [1, 2, 3], df3d)
df3d_sub21 = filter(row -> row.phase in [0, 1, 2, 3], df3d)
df3d_sub3 = filter(row -> row.phase in [4, 5, 6], df3d)


#  ==== Plots A-w curves for all phase ====
plt = @df df3d plot(
    :amp,
    :freq,
    group=:phase,
    palette=:RdYlGn_7,
    m=(2, 4),
    legend_title="phase",
    legend_position=:outertopright,
    # xlabel = "Switching Amplitude",
    # ylabel = "Switching Frequency"
    dpi=500,
    # bg = RGB(0.2, 0.2, 0.5)
)
xlabel!(plt, "Switching Amplitude")
ylabel!(plt, "Switching Frequency")
savefig(plt, save_path * "switching_dynamics_A_w_phases.png")


## === plot A-w curves for separate phase ====
#
# scheme = cgrad(:RdYlGn_11, 11, categorical = true)
plt1 = @df df3d_sub1 plot(
    :amp,
    :freq,
    group=:phase,
    # palette = scheme,
    c=["#007cd4"],
    m=(2, 4),
    legend_title="phase",
    legend_position=:topright,
    # xlims = [0, 30],
    xlabel="Switching Amplitude",
    ylabel="Switching Frequency",
    foreground_color_legend=nothing,
    dpi=500,
    # bg = RGB(0.2, 0.2, 0.5)
)
# savefig(plt1, "./figures/A-w_phase0.png")
savefig(plt1, save_path * "./A-w_phase0.png")

plt2 = @df df3d_sub2 plot(
    :amp,
    :freq,
    group=:phase,
    # palette = scheme,
    color=["#ff6232" "#ffaa4f" "#ffdf7d"],
    m=(2, 4),
    legend_title="phase",
    legend_position=:topright,
    # xlims = [0, 30],
    xlabel="Switching Amplitude",
    ylabel="Switching Frequency",
    foreground_color_legend=nothing,
    dpi=500,
    # bg = RGB(0.2, 0.2, 0.5)
)
# savefig(plt2, "./figures/A-w_phase1-3.png")
savefig(plt2, save_path * "./A-w_phase1-3.png")
#007cd4

plt21 = @df df3d_sub21 plot(
    :amp,
    :freq,
    group=:phase,
    # palette = scheme,
    color=["#007cd4" "#6541b2" "#ba3a4c" "#b38400"],
    m=(2, 4),
    legend_title="phase",
    legend_position=:topright,
    # xlims = [0, 30],
    xlabel="Switching Amplitude",
    ylabel="Switching Frequency",
    foreground_color_legend=nothing,
    dpi=500,
    # bg = RGB(0.2, 0.2, 0.5)
)
# savefig(plt21, "./figures/A-w_phase0-3.png")
savefig(plt21, save_path * "./A-w_phase0-3.png")

plt3 = @df df3d_sub3 plot(
    :amp,
    :freq,
    group=:phase,
    # palette = scheme,
    color=["#97db57" "#3ec057" "#006a31"],
    m=(2, 4),
    legendtitle="phase",
    legend_position=:topright,
    # xlims = [0, 30],
    xlabel="Switching Amplitude",
    ylabel="Switching Frequency",
    foreground_color_legend=nothing,
    dpi=500,
    # bg = RGB(0.2, 0.2, 0.5)
)
# savefig(plt3, "./figures/A-w_phase4-6.png")
savefig(plt3, save_path * "./A-w_phase4-6.png")
# plt_combo = plot(plt1, plt2, plt3, layout = (1, 3), size = (1800, 600))
# savefig(plt_combo, "phase_shift_combo.png")



## * ===## =============== A-w-ϕ curve if ϕ is not controllable================
switch_amplitude = []
switch_frequency = []
T_init = 1e-10
tspan = (0.0, 600.0)
ΔT = 100.0
@showprogress for freq in 0.0:0.1:1#exp10.(-4:0.05:1)
    for amplitude = 0:10:2000#15:0.1:20#14:0.01:18#10:0.05:30
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
            ts, cb = make_cb([T_init, T_init + ΔT], 13, amplitude)
            prob1 = ODEProblem(model, u0map, tspan, pmap)
            sol1 = solve(prob1, Rosenbrock23(), callback=cb, tstops=ts)

            # plt1 = plot(sol1, vars = [5, 6, 9, 10], lw = 1.5, title = "Amp: $amplitude, Frequency: $freq")
            # display(plt1)

            check = check_switching(sol1, ts, tspan)
            append!(switch_set, check)
            # if check == -1 # if there exist a phase
            #     append!(switch_amplitude, amplitude)
            #     append!(switch_frequency, freq)
            #     # append!(switch_phase, phase)
            # end
            if check_switching(sol1, ts, tspan) == -1
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



#
plt = plot(switch_amplitude, switch_frequency, seriestype=:scatter, #yaxis=:log10,
    label="switching events",
    # title = "Frequency vs Amplitude at the switching",
    title="ΔT = $ΔT",
    xlabel="Switch Amplitude", ylabel="Switch Frequency", dpi=500)
add_Boundary_event(switch_amplitude, switch_frequency, plt)
# savefig(plt, save_path * "switching_dynamics_A_w_independent_phases.png")


## ============ find switching time different between Dll1 and Dll4 ligand for all datasets 🔴have not finished

switch_amplitude = []
switch_frequency = []
switch_phase = []
switch_time = []

T_init = 100.0
ΔT = 200.0
tspan = [0.0, 400.0]

# for db_idx = 1:nrow(db) # iterate over each dataset row
#     @show db_idx
for amplitude = 300#:10:300
    for freq in 0.0:0.1:0.3
        single_sol = single_solve(; db_idx=db_idx, freq=freq, phase=0, amplitude=amplitude, T_init=T_init, ΔT=ΔT, tspan=tspan)
        # plt = plot(single_sol, vars = [5, 8, 6, 9], lw = 2)
        # display(plt)
        check = check_switching(single_sol)
        @show check
        t_switching = switching_time(; pulse_period=T_init:0.1:T_init+ΔT, idx=[6, 9], return_plot=false, plt=plt)
        # append!(switch_set, check)
        if check == -1
            # calculate switching time for Dll1 and Dll4 at the same time
            # note that freq = 0 is the Dll4 case.
            # t_switching  = switching_time(;pulse_period = T_init:0.1:T_init + ΔT , idx = [6,9], return_plot = false, plt = plt)
            # @show t_switching
            append!(switch_amplitude, amplitude)
            append!(switch_frequency, freq)
            append!(switch_time, t_switching)
            # append!(switch_phase, phase)

        end
    end
end

# end




## ===========Calculate the A ω ϕ st relation for a multi gene. specify the db_idx number for gene set.
# for single gene case | 49, 69
db_idx = 592
df4d = A_ω_ϕ_st_relation(amplitude_range=300, freq_range=0:0.1:2, db_idx=db_idx, phase_sample_size=15)
fixed_amp = 300
df4d_amp_300 = filter(row -> row.amp == fixed_amp, df4d)

plt_49 = plot_min_ST_ω(df4d_amp_300, plot_ontop=false)
plt_592 = plot_min_ST_ω(df4d_amp_300; figure=plt_49, plot_ontop=true, fixed_amp=fixed_amp)

##

# for multi gene case ---(ϕ indenpent)
function A_ω_ϕ_st_relation_multi_gene(; model=model, amplitude_range=0:50:300, freq_range=0:0.01:2, gene_set_id=[49, 592], phase_sample_size=6, prc2="NA", mute_parameter_disp=false)
    gene_set_A_ω_ϕ_st_df = []
    for gene_id in gene_set_id
        gene_i_df = A_ω_ϕ_st_relation(amplitude_range=amplitude_range, freq_range=freq_range, db_idx=gene_id, phase_sample_size=phase_sample_size, prc2=prc2, mute_parameter_disp=mute_parameter_disp)
        if isempty(gene_i_df)
            println("No data for gene_id = ", gene_id)
            break
        else
            @show gene_i_df
        end
        ST_ω_df = ST_ω(gene_i_df)
        ST_ω_df[!, :gene_id] .= gene_id
        @show ST_ω_df
        push!(gene_set_A_ω_ϕ_st_df, ST_ω_df)
    end
    @show gene_set_A_ω_ϕ_st_df
    return vcat(gene_set_A_ω_ϕ_st_df...)
end

df_stack = A_ω_ϕ_st_relation_multi_gene(amplitude_range=0:50:300, freq_range=0:0.01:2, gene_set_id=[49, 592], phase_sample_size=6, prc2=0.3, mute_parameter_disp=true)

# show amplitude denpendent ST- ω relation for two genes.
function plot_ST_ω_by_amplitude(; data=df_stack, layout=(3, 2), size=(1200, 800), legendfontsize=8, legendtitlefontsize=8, display_plt=false, prc2="NA")
    fig_set = []
    for amplitude in unique(data.amp)
        fixed_A_2genes = filter(row -> row.amp == amplitude, data)
        plt = @df fixed_A_2genes plot(
            :freq,
            :stime_minimum,
            group=:gene_id,
            palette=:tab10,
            m=(0.8, 1.5),
            # legend_title="Amplitude",
            legend_position=:outertopright,
            legendfontsize=legendfontsize,
            legendtitlefontsize=legendtitlefontsize,
            ylabel="Switching Time",
            xlabel="Switching Frequency",
            dpi=500,
            legend_title="Amplitude = $amplitude, \n PRC2  = $prc2"
        )
        push!(fig_set, plt)
    end
    plt = plot(fig_set..., layout=layout, size=size)
    if display_plt == true
        display(plt)
    end
    # display(plt) 
    return fig_set, plt
end

plot_ST_ω_by_amplitude()


## now test for how prc2 changes will affect 2 genes ST_ω relation in a amplitude denpendent manner
# ---- generate dataframe for 2 genes (49, 592) as prc2 rate increases ------
function Gen_df_stack_prc2_increase(; model=model, amplitude_range=0:50:300, freq_range=0:0.01:2, gene_set_id=[49, 592], phase_sample_size=10, prc2_range=0:0.1:1, mute_parameter_disp=true)
    df_stack_prc2_set = []
    @showprogress for prc2 in prc2_range
        df_stack_prc2_i = A_ω_ϕ_st_relation_multi_gene(amplitude_range=amplitude_range, freq_range=freq_range, gene_set_id=gene_set_id, phase_sample_size=phase_sample_size, prc2=prc2, mute_parameter_disp=mute_parameter_disp)
        if isempty(df_stack_prc2_i)
            println("No data for prc2 = ", prc2)
            continue
        else
            df_stack_prc2_i[!, :prc2] .= prc2
            @show df_stack_prc2_i
        end
        push!(df_stack_prc2_set, df_stack_prc2_i)
    end
    return vcat(df_stack_prc2_set...)
end



# ======= generate data for two genes with various prc2=======
df_stack_prc2_set = Gen_df_stack_prc2_increase(; model=model_pulsatile, amplitude_range=400, freq_range=0:0.1:2, gene_set_id=[49, 592], phase_sample_size=15, prc2_range=0.3:0.001:0.4, mute_parameter_disp=true)
CSV.write("df_stack_prc2_set_amplitude_range=400, freq_range=0:0.1:2, gene_set_id=[49, 592], phase_sample_size=15, prc2_range=0.3:0.001:0.4.csv", df_stack_prc2_set)
df_stack_prc2_set = CSV.File("df_stack_prc2_set_amplitude_range=300, freq_range=0:0.1:2, gene_set_id=[49, 592], phase_sample_size=15, prc2_range=0.3:0.001:0.4.csv") |> DataFrame
show(df_stack_prc2_set, allrows=true)
unique(df_stack_prc2_set.gene_id)
df_stack_prc2_set_new_gene_name = df_stack_prc2_set
# filter!(row -> row.amp == 300, df_stack_prc2_set_new_gene_name)

# generate a more detailed data for finner prc2 values.
# df_stack_prc2_finner_set = Gen_df_stack_prc2_increase(; amplitude_range=50:50:300, freq_range=0:0.05:2, gene_set_id=[49, 592], phase_sample_size=15, prc2_range=0.3:0.1:0.6, mute_parameter_disp=true)
# df_stack_prc2_finner_set
# CSV.write("df_stack_prc2_finner_set_amplitude_range=50:50:300, freq_range=0:0.05:2, gene_set_id=[49, 592], phase_sample_size=15, prc2_range=0.3:0.01:0.6.csv", df_stack_prc2_finner_set)


# =============== save and laod thie dataframe to a .csv file
using CSV
# CSV.write("df_stack_prc2_set_amplitude_range=0:50:400, freq_range=0:0.01:2, gene_set_id=[49, 592], phase_sample_size=10, prc2_range=0:0.1:1.csv", df_stack_prc2_set)
df_stack_prc2_set = CSV.File("df_stack_prc2_set_amplitude_range=0:50:400, freq_range=0:0.01:2, gene_set_id=[49, 592], phase_sample_size=10, prc2_range=0:0.1:1.csv") |> DataFrame

# ============ the loaded dataframe has gene_id with 49 and 592, I just want to map this column to "gene_1" and "gene_2"
df_stack_prc2_set_new_gene_name = df_stack_prc2_set
idx_gne_dict = Dict(49 => "Gene 1", 592 => "Gene 2")
df_stack_prc2_set_new_gene_name = transform(df_stack_prc2_set_new_gene_name, :gene_id => ByRow(x -> idx_gne_dict[x]) => :gene_id)
unique(df_stack_prc2_set_new_gene_name.gene_id)
gdf_stack_prc2 = groupby(df_stack_prc2_set_new_gene_name, :prc2)

plot_ST_ω_by_amplitude(; data=gdf_stack_prc2[3], layout=(3, 3), size=(2200, 1600), display_plt=true)




# ----compare ST ω relation only at fixed amplitude
# check prc2 values 
function select_df_prc2_fixed_amp(; data=gdf_stack_prc2, prc2_key=3, amp_select=300)
    df_prc2 = data[keys(data)[prc2_key]]
    df_prc2_fixed_amp = filter(row -> row.amp == amp_select, df_prc2)
    println("returning prc2 = ", keys(data)[prc2_key], " and amp = ", amp_select)
    return df_prc2_fixed_amp
end





keys(gdf_stack_prc2)[9]
fixed_amp = 300
df_stack_prc2_fixed_amp = select_df_prc2_fixed_amp(prc2_key=3, amp_select=fixed_amp)
fig_set, plt = plot_ST_ω_by_amplitude(data=df_stack_prc2_fixed_amp, layout=(1, 1), size=(900, 600))
plt


## ====== fixed amplitude  change prc2 for two genes =================
plt_set = []
fixed_amp = 300
for prc2_idx in 1:length(keys(gdf_stack_prc2))
    @show prc2_idx
    df_stack_prc2_fixed_amp = select_df_prc2_fixed_amp(prc2_key=prc2_idx, amp_select=fixed_amp)
    fig_set, plt = plot_ST_ω_by_amplitude(data=df_stack_prc2_fixed_amp,
        layout=(1, 1), size=(900, 600),
        legendfontsize=5, legendtitlefontsize=5,
        display_plt=false, prc2=keys(gdf_stack_prc2)[prc2_idx][1])
    # title!(plt, "PRC2 = $(keys(gdf_stack_prc2)[prc2_idx][1])", titlefont=font(8))
    xlims!(plt, (0, 2))
    ylims!(plt, (100, 150))
    push!(plt_set, plt)
end

plt_ST_ω_2gene_prc2_increase_fixed_amp = plot(plt_set[1:4]..., layout=(2, 2), size=(1600, 1000), legendfontsize=12, legendtitlefontsize=12, m=(3, 0.6), legend=:topright)
# savefig(plt_ST_ω_2gene_prc2_increase_fixed_amp, "plt_ST_ω_2gene_prc2_increase_fixed_amp.png")




# ========== make an animation for prc2 0.3 ~ 0.4 at amp = 300 for two genes (49,592) ==================
plt_set = []
fixed_amp = 300
anim_prc2 = @animate for prc2_idx in 1:length(keys(gdf_stack_prc2))
    @show prc2_idx
    df_stack_prc2_fixed_amp = select_df_prc2_fixed_amp(prc2_key=prc2_idx, amp_select=fixed_amp)
    fig_set, plt = plot_ST_ω_by_amplitude(data=df_stack_prc2_fixed_amp,
        layout=(1, 1), size=(900, 600),
        legendfontsize=5, legendtitlefontsize=5,
        display_plt=false, prc2=keys(gdf_stack_prc2)[prc2_idx][1])
    # title!(plt, "PRC2 = $(keys(gdf_stack_prc2)[prc2_idx][1])", titlefont=font(8))
    xlims!(plt, (0, 2))
    ylims!(plt, (100, 200))
    push!(plt_set, plt)
end

gif(anim_prc2, "prc2_0.3~0.4_amp_200_id(49,592).gif", fps=15)


fixed_amp = 400
gene2_sustained_df = filter(row -> row.freq == 0 && row.gene_id == "Gene 2", df_stack_prc2_set_new_gene_name)
gene1_pusatile_df = filter(row -> row.freq == 0.9 && row.gene_id == "Gene 1", df_stack_prc2_set_new_gene_name)

plt1 = @df gene1_pusatile_df plot(:prc2, :stime_minimum,
    # palette=:RdYlBu_6,
    palette=:tab10,
    # m=(0.8, 2),
    legend_title="Genes",
    legend_position=:topright,
    linewidth=2,
    xlabel="PRC2 rate",
    ylabel="Switching Time (ST)",
    dpi=500,
    foreground_color_legend=nothing,
    title="Amplitude = $fixed_amp",
    label="Gene 1 (Pusatile frequency = 0.9)",
)
plt2 = @df gene2_sustained_df plot!(plt1, :prc2, :stime_minimum,
    palette=:RdYlBu_6,
    # m=(0.8, 2),
    legend_title="Genes",
    legend_position=:topright,
    linewidth=2,
    xlabel="PRC2 rate",
    ylabel="Switching Time (ST)",
    dpi=500,
    foreground_color_legend=nothing,
    title="Amplitude = $fixed_amp",
    label="Gene 2 (Sustained)"
)
ylims!(plt2, (100, 150))
savefig(plt2, "plt2_ST_ω_2gene_prc2_increase_fixed_amp_400_freq_0.9.png")




## ========= To generate a plot of ω vs. ST for a range of PRC2 rate and a range of amplitude ==========
# ------- generate dataset
# df_stack_prc2_set_gene_49 = Gen_df_stack_prc2_increase(; model=model_pulsatile, amplitude_range=100:200:300, freq_range=0:0.1:2, gene_set_id=[49], phase_sample_size=15, prc2_range=0.1:0.05:0.6, mute_parameter_disp=true)
df_stack_prc2_set_gene_49 = Gen_df_stack_prc2_increase(; model=model_pulsatile, amplitude_range=100:200:300, freq_range=0:0.1:2, gene_set_id=[49], phase_sample_size=5, prc2_range=0.41, mute_parameter_disp=true)
# CSV.write("df_stack_prc2_set_gene_49.csv", df_stack_prc2_set_gene_49)


# ------- plot ω vs. ST for various prc2 with amplitude  = 300
df_stack_prc2_set_gene_49_fixed_amp = filter(row -> row.amp == 300, df_stack_prc2_set_gene_49)
df_stack_prc2_set_gene_49_fixed_amp
plt_id49_amp_300_prc2_set = @df df_stack_prc2_set_gene_49_fixed_amp plot(
    :freq,
    :stime_minimum,
    group=:prc2,
    palette=:RdBu_9,
    m=(0.8, 1.5),
    # legend_title="Amplitude",
    legend_position=:outertopright,
    legendfontsize=5,
    legendtitlefontsize=5,
    ylabel="Switching Time",
    xlabel="Switching Frequency",
    dpi=500,
    legend_title="Amplitude = 300, \n PRC2 rate"
)
# savefig(plt_id49_amp_300_prc2_set, "plt_id49_amp_300_prc2_set.png")



# ------- plot ω vs. ST for various prc2 with amplitude  = 100
df_stack_prc2_set_gene_49_fixed_amp = filter(row -> row.amp == 100, df_stack_prc2_set_gene_49)
df_stack_prc2_set_gene_49_fixed_amp
plt_id49_amp_100_prc2_set = @df df_stack_prc2_set_gene_49_fixed_amp plot(
    :freq,
    :stime_minimum,
    group=:prc2,
    palette=:RdBu_9,
    m=(0.8, 1.5),
    # legend_title="Amplitude",
    legend_position=:outertopright,
    legendfontsize=5,
    legendtitlefontsize=5,
    ylabel="Switching Time",
    xlabel="Switching Frequency",
    dpi=500,
    legend_title="Amplitude = 100, \n PRC2 rate"
)
# savefig(plt_id49_amp_100_prc2_set, "plt_id49_amp_100_prc2_set.png")


# ----- compared with paper figure 6 with default prc2 rate 0.41--------------------------------
df4d = A_ω_ϕ_st_relation(; model=model_pulsatile,
    amplitude_range=100:200:300, freq_range=0:0.02:2,
    ΔT=100, db_idx=49)
ST_ω_df = ST_ω(df4d)
df4d_amp_300 = filter(row -> row.amp == 300, ST_ω_df)
@df df4d_amp_300 plot(
    :freq,
    :stime_minimum,
    # group=:phase,
    palette=:RdBu_6,
    m=(0.8, 1.5),
    # legend_title="Amplitude",
    legend_position=:outertopright,
    legendfontsize=5,
    legendtitlefontsize=5,
    ylabel="Switching Time",
    xlabel="Switching Frequency",
    dpi=500,
    legend_title="Amplitude =300, \n PRC2 rate = 0.41"
)


