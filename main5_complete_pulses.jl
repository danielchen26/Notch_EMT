# This file is a modified version of main4 file which fixed the problem of not having a complete pulse during a certain duration of time of input signal. 



## 🥗====================== Loading packages and data library==========================
include("./main5_pkgs.jl")
include("./main5_functions.jl")
using Revise
db, p_names, initi_names, parameter_set, initial_condition, rand_idx_set, db_idx = loading_database()




## 🥗============================= Import models and signal function =================================
signal_type = "bump"
model_bump = Import_model(; type=signal_type)
@show equations(model_bump)
signal_type = "pulsatile"
model_pulsatile = Import_model(; type=signal_type)
@show equations(model_pulsatile)
rn_latex, ode_latex = ode_latex_gen(model_pulsatile)
# ------------- input signal function --------------------------------
osci_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)))



## 🥗 ============================= a single run for a constant Dll4 signal ==============================
# ======= Dll4 signal parameters
signal = Signal(
    db_idx=592,
    tspan=(0.0, 350.0),
    freq=0.0,
    amplitude=220.0,
    phase=0.0, # if constant DLL4 signal is given, then phase has to be zero.
    T_init=100.0,
    ΔT_Dll4_ref=100.0,
    ΔT=100.0)
# reset_signal(; signal=signal::Signal, amplitude=signal.amplitude, freq = signal.freq)
# ======= solve
u0map, pmap, p, tspan, ts, cb, sol = single_solve(; model=model_pulsatile, signal=signal);
# ======= plot
plt = plot(sol, idxs=[4, 5, 6, 9], lw=1.5, xlabel="Time", ylabel="Concentration", dpi=500)
tt = signal.T_init:0.01:signal.T_init+signal.ΔT
plot!(plt,
    tt,
    osci_signal.(tt, signal.amplitude, signal.freq, 0.0),
    label="Pulsatile Input", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
    # ylims = [0, 700],
    fill=(0, 0.3, :darkgreen), color="black", dpi=300)




## ====== I want to find if 🔴 changing the prc2 kinetic rate is able to result in early activation
prob_new = remake_prob(model_pulsatile, u0map, tspan, p; prc2=0.59)
@time sol = solve(prob_new, Rosenbrock23(), callback=cb, tstops=ts)
plt = plot(sol,
    idxs=[4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
    lw=1.5,
    xlabel="Time", ylabel="Concentration",
    title="PRC2 rate : 0.4",
    dpi=500)

# I need to change the signal phase every time I change the prc2 rate? / 
anim_prc2_changing(model=model_pulsatile, 0:0.1:0.8, tspan=[0, 350], u0=u0map)
reset_signal


## 🥗 =============================== a single run for a pulsatile Dll1 signal ==============================
function frequency_behaviour_fixed_amp(; id, input_freq)
    signal = Signal(
        db_idx=id,
        tspan=(0.0, 350.0),
        freq=input_freq,
        amplitude=220.0,
        phase=0.0,
        T_init=100.0,
        ΔT_Dll4_ref=150.0,
        ΔT=150.0
    )
    model = model_pulsatile
    signal_type = "pulsatile"
    plt_pulsatile = single_solve_plot(; model=model, signal=signal, type=signal_type)
    display(plt_pulsatile)
    # reset_signal(; signal = signal, amplitude=signal.amplitude, freq=signal.freq)
    # @show signal
    # # ======= solve
    # u0map, pmap, p, tspan, ts, cb, sol = single_solve(; model=model_pulsatile, signal=signal)
    # # ======= plot
    # plt = plot(sol, vars=[4, 5, 6, 9], lw=1.5, xlabel="Time", ylabel="Concentration", dpi=500)
    # tt = signal.T_init:0.01:signal.T_init+signal.ΔT# -period *.123
    # plot!(plt, tt, osci_signal.(tt, signal.amplitude, signal.freq, signal.phase),
    #     label="Pulsatile Input", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
    #     # ylims = [0, 700],
    #     fill=(0, 0.3, :darkgreen), color="black", dpi=300)
    # display(plt)
end

for freq = 0.1:0.02:1
    frequency_behaviour_fixed_amp(; id=49, input_freq=freq)
end






## 🥗====== a test for gene 49 that requires multiple pulses =======
# the signal phase is adjusted.
plotly()
signal = Signal(
    db_idx=49,
    tspan=(0.0, 350.0),
    freq=0.1,
    amplitude=57.5,
    phase=0,
    T_init=100,
    ΔT_Dll4_ref=150.0,
    ΔT=150.0
)
reset_signal(; signal, amplitude=signal.amplitude, freq=signal.freq, T_init=signal.T_init, ΔT_Dll4_ref=signal.ΔT_Dll4_ref)
signal.T_init = 70
dump(signal)
# ======= pulsatile model
model = model_pulsatile
signal_type = "pulsatile"
_, _, _, _, ts, _, single_sol = single_solve(; model=model, signal=signal);
@show check = check_switching(single_sol, [signal.T_init, signal.T_init + signal.ΔT], signal.tspan)
plt_pulsatile = single_solve_plot(; model=model, signal=signal, type=signal_type)


##
function change_t_init(T_init)
    signal = Signal(
    db_idx=49,
    tspan=(0.0, 350.0),
    freq=0.1,
    amplitude=507.5,
    phase=0,
    T_init=T_init,
    ΔT_Dll4_ref=150.0,
    ΔT=150.0
)
    reset_signal(; signal, amplitude=signal.amplitude, freq=signal.freq, T_init=signal.T_init, ΔT_Dll4_ref=signal.ΔT_Dll4_ref)
    dump(signal)
    # ======= pulsatile model
    model = model_pulsatile
    signal_type = "pulsatile"
    _, _, _, _, ts, _, single_sol = single_solve(; model=model, signal=signal);
    @show check = check_switching(single_sol, [signal.T_init, signal.T_init + signal.ΔT], signal.tspan)
    plt_pulsatile = single_solve_plot(; model=model, signal=signal, type=signal_type)
    display(plt_pulsatile)
end 


# write loop over T_init from 0 to 500
for T_init = 0:5:500
    change_t_init(T_init)
end














## 🥗================  Calculate the A ω ϕ st relation for a single gene. specify the db_idx number for a gene.========================
signal = Signal(
    db_idx=209,
    tspan=(0.0, 2550.0),
    freq=0.5,
    amplitude=65.0,
    phase=0.0,
    T_init=100.0,
    ΔT_Dll4_ref=2000.0,
    ΔT=2000.0
)
# db_idx = 600 # 🍏 paper figure 6
df4d_pulsatile = A_ω_ϕ_st_relation(; model=model_pulsatile,
    amplitude_range=0:50:300, freq_range=0.001:0.02:2,
    ΔT_Dll4_ref=100, db_idx=signal.db_idx)

plt = @df df4d_pulsatile plot(
    :freq,
    :stime,
    group=:amp,
    palette=:RdYlBu_6,
    # palette = Symbol("RdYlBu_"*"$color_catg"),
    m=(2, 4),
    legend_title="Amplitude",
    legend_position=:outertopright,
    # xlabel = "Switching Amplitude",
    # ylabel = "Switching Frequency"
    dpi=500,
    # title = "Amplitude = 100"
    # bg = RGB(0.2, 0.2, 0.5)
)
xlabel!(plt, "Driving Frequency")
ylabel!(plt, "Switching Time (ST)")


plt = plot(df4d_pulsatile.amp, df4d_pulsatile.freq, seriestype=:scatter, #yaxis=:log10,
    label="switching events", title="Frequency vs Amplitude at the switching",
    xlabel="switch_amplitude", ylabel="switch_frequency", dpi=500)
add_Boundary_event(df4d_pulsatile.amp, df4d_pulsatile.freq, plt)




















## ===## =============== A-w-ϕ curve if ϕ is not controllable================
# need to  modified the phase ---

switch_amplitude = []
switch_frequency = []
signal = Signal(
    db_idx=209,
    tspan=(0.0, 2550.0),
    freq=0.5,
    amplitude=65.0,
    phase=0.0,
    T_init=100.0,
    ΔT_Dll4_ref=2000.0,
    ΔT=2000.0
)

@showprogress for freq in 0.0:0.01:2#exp10.(-4:0.05:1)
    switch_set = []
    for amplitude = 0:1:300#15:0.1:20#14:0.01:18#10:0.05:30

        reset_signal(; signal, amplitude=amplitude, freq=freq)
        # dump(signal)
        println("\n")

        model = model_pulsatile
        signal_type = "pulsatile"
        _, _, _, _, ts, _, single_sol = single_solve(; model=model, signal=signal)
        @show check = check_switching(single_sol, [signal.T_init, signal.T_init + signal.ΔT], signal.tspan)

        check = check_switching(single_sol, ts, signal.tspan)
        append!(switch_set, check)
        if check == -1 # if there exist a phase
            append!(switch_amplitude, amplitude)
            append!(switch_frequency, freq)
        end
        if -1 in switch_set
            break
        end
    end

end


plotly()
#
plt = plot(switch_amplitude, switch_frequency, seriestype=:scatter, #yaxis=:log10,
    label="switching events",
    # title = "Frequency vs Amplitude at the switching",
    # title="ΔT = $ΔT",
    xlabel="Switch Amplitude", ylabel="Switch Frequency", dpi=500)
add_Boundary_event(switch_amplitude, switch_frequency, plt)


















## 🥗============== poisson pulses of above version =============
signal = Signal(
    db_idx=62,
    tspan=(0.0, 350.0),
    freq=0.5,
    amplitude=65.0,
    phase=0.0,
    T_init=100.0,
    ΔT_Dll4_ref=150.0,
    ΔT=150.0
)
# db_idx = 600 # 🍏 paper figure 6
df4d_poisson = A_ω_ϕ_st_relation_poisson(; amplitude_range=0:100:300, freq_range=0.001:0.2:2,
    ΔT_Dll4_ref=100, db_idx=signal.db_idx)
# CSV.write("./Data/Poisson/Poisson idx_592: freq_vs_st_amp_df: amp=0:50:300 freq = 0.001:0.1:2 ΔT = 150.csv", df4d_poisson)
plt = @df df4d_poisson plot(
    :freq,
    :stime,
    group=:amp,
    palette=:RdYlBu_6,
    # palette = Symbol("RdYlBu_"*"$color_catg"),
    m=(2, 4),
    legend_title="Amplitude",
    legend_position=:outertopright,
    # xlabel = "Switching Amplitude",
    # ylabel = "Switching Frequency"
    dpi=500,
    # title = "Amplitude = 100"
    # bg = RGB(0.2, 0.2, 0.5)
)
xlabel!(plt, "Driving Frequency")
ylabel!(plt, "Switching Time (ST)")

plt = @df df4d_poisson plot(
    :amp,
    :freq,
    # group =  :amp,
    # palette = :RdYlBu_6,
    # palette = Symbol("RdYlBu_"*"$color_catg"),
    m=(2, 4),
    # legend_title = "Amplitude",
    legend_position=:outertopright,
    # xlabel = "Switching Amplitude",
    # ylabel = "Switching Frequency"
    dpi=500,
    # title = "Amplitude = 100"
    # bg = RGB(0.2, 0.2, 0.5)
)










## 🥗=============== Poisson: Driving amp vs. driving freq ================
switch_amplitude = []
switch_frequency = []
# T_init = 100.0
# tspan = (0.0, 600.0)
# ΔT = 100.0
# db_idx = 49
# ΔT_Dll4_ref = 150.0
# signal = Signal(
#                 db_idx = 49,
#                 tspan = (0.0, 600.0),
#                 freq = 0.5, 
#                 amplitude = 65.0, 
#                 phase = 0.0, 
#                 T_init = 100.0, 
#                 ΔT_Dll4_ref = 150.0,
#                 ΔT = 100.0
#                 )
signal = Signal(
    db_idx=49,
    tspan=(0.0, 350.0),
    freq=0.01,
    amplitude=40.0,
    phase=0.0,
    T_init=100.0,
    ΔT_Dll4_ref=150.0,
    ΔT=150.0
)

# model = model_pulsatile # 🍕 for poisson model, no need to specify
@showprogress for freq in 0.001:0.2:2#exp10.(-4:0.05:1)
    switch_set = []
    for amplitude = 0:4:200#15:0.1:20#14:0.01:18#10:0.05:30
        # signal.db_idx = db_idx
        # signal.tspan = tspan
        # signal.amplitude = amplitude
        # signal.freq = freq
        # signal.phase = 3*pi/2 - signal.freq*signal.T_init # 🔴reset phase
        # signal.T_init = T_init
        # signal.ΔT_Dll4_ref = ΔT_Dll4_ref
        # period = 2pi/signal.freq
        # num_cycle= signal.ΔT_Dll4_ref/(period)
        # if floor(num_cycle) < 1.0
        #     signal.ΔT = num_cycle* period
        #     println("ΔT is :\n", signal.ΔT )
        # else
        #     signal.ΔT =  floor(num_cycle)* period  - 0.01
        #     println("ΔT is :\n", signal.ΔT )
        # end 

        @show freq, amplitude
        # ---- 🍕 regular pulses -----
        # reset_signal(signal = signal, freq = freq, amplitude = amplitude)
        # @show signal
        # _, _, _, _, ts, _, single_sol = single_solve(; model=model, signal = signal::Signal, prc2="NA")




        # ---- 🍕 poisson pulses -----
        signal.amplitude = amplitude
        signal.freq = freq
        @show signal
        u0map, pmap, p, signal.tspan, ts, cb, single_sol = single_solve_poisson(; signal=signal::Signal, prc2="NA", mute_parameter_disp=true)
        # plt1 = plot(sol1, vars = [5, 6, 9, 10], lw = 1.5, title = "Amp: $amplitude, Frequency: $freq")
        # display(plt1)

        # check = check_switching(single_sol, ts, signal.tspan)
        check = check_switching(single_sol, ts, signal.tspan)
        append!(switch_set, check)
        # @show signal
        @show check
        if check == -1 # if there exist a phase
            append!(switch_amplitude, amplitude)
            append!(switch_frequency, freq)
            # append!(switch_phase, phase)
        end
        # if check_switching(single_sol, ts, signal.tspan) == -1
        #     break
        # end

        # @show switch_set
        # if -1 in switch_set
        #     append!(switch_amplitude, amplitude)
        #     append!(switch_frequency, freq)
        # end
        if -1 in switch_set
            break
        end

    end
end
freq_vs_amp_df = DataFrame(freq=switch_frequency, amp=switch_amplitude)
CSV.write("./Data/Poisson/Poisson: freq_vs_amp_df: freq = 0.001:0.2:2 amp= 0:4:200 ΔT = 150 Ins_2.csv", freq_vs_amp_df)
plt = plot(switch_amplitude, switch_frequency, seriestype=:scatter, #yaxis=:log10,
    label="switching events",
    # title = "Frequency vs Amplitude at the switching",
    legend_position=:topright,
    title="ΔT = $(signal.ΔT)",
    xlabel="Driving Amplitude", ylabel="Driving Frequency", dpi=500)
add_Boundary_event(switch_amplitude, switch_frequency, plt)













## 🥗=========================================  Poisson distributed pulses case  =====================================================
# db_idx = 592
# tspan = (0.0, 350.0)
# T_init = 100.0 
# ΔT = 150.0
# freq = 0.2
# amplitude = 365.0

# # ====== each time loading poisson model will generate a sequence of poisson pulses within [T_init, T_init + ΔT]
# model, single_instance_sequence = Load_poisson_model(freq,T_init,ΔT) 
# p = vcat([collect(parameter_set[db_idx, :]),  0.0]...)
# pmap = parameters(model) .=> p
# u0 = collect(initial_condition[db_idx, :])
# u0map = species(model) .=> u0
# ts, cb = make_cb([T_init, T_init + ΔT], 12, amplitude)

# prob1 = ODEProblem(model, u0map, tspan, pmap)
# sol = solve(prob1, Rosenbrock23(), callback=cb, tstops=ts);





signal = Signal(
    db_idx=49,
    tspan=(0.0, 350.0),
    freq=0.01,
    amplitude=40.0,
    phase=0.0,
    T_init=100.0,
    ΔT_Dll4_ref=150.0,
    ΔT=150.0
)
u0map, pmap, p, signal.tspan, ts, cb, sol = single_solve_poisson(; signal=signal::Signal, prc2="NA", mute_parameter_disp=true)
check = check_switching(sol, ts, signal.tspan)
plot(sol,
    vars=[4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
    lw=2,
    xlabel="Time", ylabel="Concentration",
    foreground_color_legend=nothing,
    dpi=500)




## ========== pulses shape  ========
# regular square pulse shape
freq = 1.5
T_init = 100
ΔT = 100
single_instance_sequence = generate_poisson_pulses(freq, T_init, T_init + ΔT)
FR(t) = compose_pulses(t, single_instance_sequence)
plot(0:0.01:300, FR.(0:0.01:300))



## === hill function round pulse =================================

function interval_rd(t, t_middle, shift, n, k)
    hill(t + t_middle - shift, n, k) - hill(t - t_middle - shift, n, k)
end


function hill(t, n, k)
    return (t)^n / ((t)^n + k^n)
end

function compose_pulses_rd(t, intervals, shift, n, k)
    result = 0.0
    for i in 1:length(intervals)
        t_middle = (intervals[i][1] + intervals[i][2]) / 2
        # width =  K^(1/n) - (K/(2^(1/n)-1))^(1/n)
        result += interval_rd(t, t_middle, shift, n, k)
    end
    return result
end

tt = collect(0:0.01:13) .+ 3
# interval_rd.(tt, 5, 10, 1)
f(t) = compose_pulses_rd(t, [[4, 6], [8, 10]], 1, 100, 1)
plot(tt, f.(tt))
xticks!(1:0.5:13)

tt = collect(0:0.01:13) .+ 3
plot(tt, hill.(tt, 100, 1))
xticks!(0:1:13)









## =============================================================================
## ======= Two plots Dll4 vs Dll1 for gene id:49 ====
signal1 = Signal(
    db_idx=49,
    tspan=(0.0, 350.0),
    freq=0.25,
    amplitude=40.0,
    phase=0.0,
    T_init=100.0,
    ΔT_Dll4_ref=150.0,
    ΔT=150.0
)
prc2 = 0.41



# 🍕 
reset_signal(; signal=signal1, amplitude=signal1.amplitude, freq=signal1.freq)
dump(signal1)


plt_gene1_Dll4, plt_gene1_Dll1, plt_2genes_compare_id_49 =
    Two_Genes_TS_by_Prc2(;
        model=model_pulsatile,
        signal1=signal1, signal2=signal1, prc2=prc2,
        title_on=true, legend_title_on=false,
        vars_to_show=[5, 6, 9], #tspan_rt = 2, # 
        type="pulsatile") #🔴 specify the type
plt_gene1_Dll4
plt_gene1_Dll1
plt_2genes_compare_id_49
# savefig(plt_2genes_compare_id_49,"./figures/APS_id_49_Dll4_Dll1_compare_bump.png")


p = vcat([collect(parameter_set[signal1.db_idx, :]), signal1.freq, 0.0, signal1.phase]...)
@show p
u0 = collect(initial_condition[signal1.db_idx, :])
u0map = species(model) .=> u0
@show u0map
@show remake_prob(model, u0map, signal1.tspan, p; prc2=0.4, mute_parameter_disp=false)



remake_single_solve(; model=model, signal=signal1, prc2=prc2, mute_parameter_disp=false)