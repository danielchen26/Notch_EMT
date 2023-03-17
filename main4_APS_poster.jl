## ====================== Loading packages and data library==========================
using DifferentialEquations
using Plots
using Catalyst
using Catalyst: parameters
using ProgressMeter
using DataFrames, CSV, Random, Distributions
using StatsPlots
using Latexify, Measures
include("./Functions.jl")
db, p_names, initi_names, parameter_set, initial_condition, rand_idx_set, db_idx = loading_database()



## ====================== Load model =================================================
signal_type = "bump"
model_bump = Import_model(; type = signal_type)
@show equations(model_bump)
signal_type = "pulsatile"
model_pulsatile = Import_model(; type = signal_type)
@show equations(model_pulsatile)



## ======= Two plots Dll4 vs Dll1 for gene id:49 ====
single_gene_id = 49
id2_freq = 0.25
phase2 = 0
amplitude1 = 65
amplitude2 = 65
T_init = 20
ΔT = 100
prc2 = 0.41

plt_gene1_Dll4, plt_gene1_Dll1, plt_2genes_compare_id_49 =
    Two_Genes_TS_by_Prc2(;
        model = model_bump,
        # model = model_pulsatile,
        id1 = single_gene_id, id2 = single_gene_id,
        id2_freq = id2_freq, phase2 = phase2, amplitude1 = amplitude1,
        amplitude2 = amplitude2, prc2 = prc2,
        T_init = T_init, ΔT = ΔT, title_on = true, legend_title_on = true,
        vars_to_show = [5, 6, 9], #tspan_rt = 2, # 
        type = "bump") #🔴 specify the type
plt_gene1_Dll4
plt_gene1_Dll1
plt_2genes_compare_id_49



## ======= Two plots 🔴 animation Dll4 vs Dll1 for 1 gene ====
single_gene_id = 49
id2_freq = 0.1
phase2 = 3
amplitude1 = 60
amplitude2 = 100
prc2 = 0.41
Two_genes_prc2_TS_animation(; prc2_range = 0:0.02:1, model = model_pulsatile,
    db_idx1 = single_gene_id, db_idx2 = single_gene_id,
    id2_freq = 0.3, phase2 = phase2,
    amplitude1 = amplitude1, amplitude2 = amplitude1, type = "pulsatile"
)

## visualization for the two genes TS by changing prc2 rate

_, _, plt_2genes_compare = Two_Genes_TS_by_Prc2(; model = model_bump,
    # id1 = 49, id2 = 593,
    id1 = 593, id2 = 49,
    id2_freq = 0.1, prc2 = 0.192,
    phase2 = 0,
    # T_init = 0,
    amplitude1 = 120, amplitude2 = 120,
    type = "bump")
plt_2genes_compare
# savefig(plt_2genes_compare, "./figures/592_49_prc2=0.862.png")

## === two gene bump model changing prc2 rate ==== 
Two_genes_prc2_TS_animation(; model = model_bump, prc2_range = 0:0.02:1,
    db_idx1 = 593, db_idx2 = 49,
    amplitude1 = 100, amplitude2 = 100,
    type = "bump")
Two_genes_prc2_TS_animation(; model = model_bump, prc2_range = 0:0.02:1,
    db_idx1 = 49, db_idx2 = 593,
    amplitude1 = 90, amplitude2 = 90,
    type = "bump")

## === two gene pulsatile model changing prc2 rate ==== 
Two_genes_prc2_TS_animation(; model = model_pulsatile, prc2_range = 0:0.02:1,
    db_idx1 = 593, db_idx2 = 49, phase2 = 0,
    amplitude1 = 130, amplitude2 = 130,
    type = "pulsatile")
Two_genes_prc2_TS_animation(; model = model_pulsatile, prc2_range = 0:0.02:1,
    db_idx1 = 49, db_idx2 = 593, phase2 = 1,
    amplitude1 = 130, amplitude2 = 130,
    type = "pulsatile")




## =========  Calculate the
db_idx = 49
df4d = A_ω_ϕ_st_relation(; model = model_pulsatile,
    amplitude_range = 0:50:300, freq_range = 0:0.02:2,
    ΔT = 100, db_idx = db_idx)

df4d_phase_0, _, _, df4d_phase_3, _, _, df4d_phase_6, df4d_amp_select = df4d_sub_gen(df4d, amplitude_select = 300)

plt_fix_phase_freq_vs_ST = df4d_fix_phase_freq_vs_ST(df4d_phase_6; ΔT = 100, save = true)
plt_fix_amplitude_freq_vs_ST = df4d_fix_amp_freq_vs_ST(df4d_amp_select; ΔT = 100, save = true)


## =======
switch_amplitude = []
switch_frequency = []
switch_phase = []
tspan = [0.0, 150.0]
# model = model_bump
model = model_pulsatile # model_bump
model = model_bump # model_bump
ΔT = 100
@showprogress for phase = 0:2π
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
            sol1 = solve(prob1, Rosenbrock23(), callback = cb, tstops = ts)

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

plt = plot(switch_amplitude, switch_frequency, seriestype = :scatter, #yaxis=:log10,
    label = "switching events", title = "Frequency vs Amplitude at the switching",
    xlabel = "switch_amplitude", ylabel = "switch_frequency", dpi = 500)
add_Boundary_event(switch_amplitude, switch_frequency, plt)

df3d = DataFrame(freq = switch_frequency, phase = switch_phase, amp = switch_amplitude)
df3d_sub1 = filter(row -> row.phase in [0], df3d)
df3d_sub2 = filter(row -> row.phase in [1, 2, 3], df3d)
df3d_sub21 = filter(row -> row.phase in [0, 1, 2, 3], df3d)
df3d_sub3 = filter(row -> row.phase in [4, 5, 6], df3d)
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