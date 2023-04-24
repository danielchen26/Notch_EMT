
# module main5_functions
# export Signal, loading_database, Import_model, ode_latex_gen
## ====================== Defining the model and parameters ==========================
# module model_parameters
#     export Signal
#     Base.@kwdef mutable struct Signal
#         db_idx::Int64
#         tspan::Tuple{Float64, Float64}
#         freq::Float64
#         phase::Float64
#         amplitude::Float64
#         T_init::Float64
#         ΔT_Dll4_ref ::Float64 # this is the reference ΔT for DLL4, for each frequency in Dll1, the ΔT will have some adjustment to ensure the at least two complete pulses are generated.
#         ΔT::Float64
#     end
# end

Base.@kwdef mutable struct Signal
    db_idx::Int64
    tspan::Tuple{Float64,Float64}
    freq::Float64
    phase::Float64
    amplitude::Float64
    T_init::Float64
    ΔT_Dll4_ref::Float64
    ΔT::Float64
end

Base.@kwdef mutable struct Poisson_Signal
    db_idx::Int64
    tspan::Tuple{Float64,Float64}
    freq::Float64
    amplitude::Float64
    T_init::Float64
    ΔT_Dll4_ref::Float64
    ΔT::Float64
end






## ================ import database =================
function loading_database(; data_path="../../../Notch_EMT_data/Notch_params_complete.csv")
    db = CSV.File(data_path) |> DataFrame
    p_names = names(db)[1:11]
    initi_names = names(db)[13:end]

    # df = db[shuffle(1:nrow(db))[1:10], :]
    # parameter_set = df[:, p_names]
    # initial_condition = df[:, initi_names]
    ##  ================================ With random parameters generate switching dynamics ============================================
    # db_idx = 210 is very sensitive to low amplitude = 0.9
    # df = db[shuffle(1:nrow(db))[1:10], :]
    parameter_set = db[:, p_names]
    initial_condition = db[:, initi_names]
    rand_idx_set = shuffle(1:nrow(db))
    @show db_idx = rand(rand_idx_set)
    @show parameter_set[db_idx, :]
    # db_idx = 1961
    return db, p_names, initi_names, parameter_set, initial_condition, rand_idx_set, db_idx
end

##  ============ Import models ==========

function Import_model(; type="pulsatile")
    if type == "pulsatile"
        model_pulsatile = @reaction_network begin
            # @parameter k0 k1 k2 d m p k pp kk δ α1 w A ϕ # put A at last as the control 13th variable
            # @species M(t) R(t) MR(t)
            (A * (1 + sign(cos(w * t + ϕ))), 1.0), R ↔ NR               # NICD binds RBPJ
            # (A * (abs(cos(w * t + ϕ))), 1.0), R ↔ NR               # NICD binds RBPJ 🍏🔴
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
        model = model_pulsatile
    elseif type == "bump"
        model_bump = @reaction_network begin
            # @parameters k0 k1 k2 d m p k pp kk δ α1 w A ϕ # put A at last as the control 13th variable
            # (A * (1 + sign(cos(w * t + ϕ))), 1.0), R ↔ NR               # NICD binds RBPJ
            (A * (abs(cos(w * t + ϕ))), 1.0), R ↔ NR               # NICD binds RBPJ 🍏🔴
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
        end k0 k1 k2 d m p k pp kk δ α # put A at last as the control 13th variable
        model = model_bump
    end
    @show species(model)
    @show parameters(model)
    @show speciesmap(model)
    @show paramsmap(model)
    return model
end


function ode_latex_gen(model::ReactionSystem)
    rn_latex = latexify(model)
    ODE_equations = convert(ODESystem, model)
    ode_latex = latexify(ODE_equations)
    return rn_latex, ode_latex
end


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
    cb = DiscreteCallback(condition, affect!, save_positions=(true, true))
    return ts, cb
end

"""
test if the steady state switched
	-1 : switched
	 1 : not swithced
"""
function check_switching(sol, ts, tspan)
    H4_i, H27_i = sol(ts[1])[[6, 9]]
    H4_f, H27_f = sol(tspan[2])[[6, 9]]
    init = H4_i / H27_i < 1 ? 1 : -1
    after = H4_f / H27_f < 1 ? 1 : -1
    return init * after
end

function add_Boundary_event(switch_amplitude, switch_frequency, plt)
    df = DataFrame(switch_amplitude=switch_amplitude, switch_frequency=switch_frequency)
    # CSV.write("switching_freq_amplitude.csv",df)

    gp = groupby(df, :switch_frequency)
    boundary_amplitude = []
    for i in gp
        append!(boundary_amplitude, minimum(i.switch_amplitude))
    end

    critical_freq = [keys(gp)[i].switch_frequency for i = 1:gp.ngroups]
    boundary_amplitude
    plt_boundary = plot!(plt, boundary_amplitude, critical_freq, label="switching boundary", lw=3)
end

# save A-w data
function save_A_w_data(switch_amplitude, switch_frequency, path)
    df = DataFrame(switch_amplitude=switch_amplitude, switch_frequency=switch_frequency)
    CSV.write(path * "switch_amplitude.csv", df)
end


# save initial condition, parameter set data
function save_init_param(parameter_set, initial_condition, id, path)
    CSV.write(path * "initial_condition.csv", DataFrame(initial_condition[id, :]))
    CSV.write(path * "parameter_set.csv", DataFrame(parameter_set[id, :]))
end


# single solve, given a instance in database specified by db_idx.
function single_solve(; model=model, signal=signal::Signal, prc2="NA", mute_parameter_disp=false) #🍏🔴added prc2 changing option
    p = vcat([collect(parameter_set[signal.db_idx, :]), signal.freq, 0.0, signal.phase]...)
    pmap = parameters(model) .=> p
    u0 = collect(initial_condition[signal.db_idx, :])
    u0map = species(model) .=> u0
    ts, cb = make_cb([signal.T_init, signal.T_init + signal.ΔT], 13, signal.amplitude)
    prob1 = ODEProblem(model, u0map, signal.tspan, pmap)
    if prc2 != "NA"
        prob1 = remake_prob(model, u0map, signal.tspan, p; prc2=prc2, mute_parameter_disp=mute_parameter_disp)
    end
    sol = solve(prob1, Rosenbrock23(), callback=cb, tstops=ts)
    # sol = solve(prob1, Rosenbrock23(), callback=cb, tstops=ts, abstol=1e-10, reltol=1e-8)
    # sol = solve(prob1, TRBDF2(), callback=cb, tstops=ts);
    return u0map, pmap, p, signal.tspan, ts, cb, sol
end


# write a function to remake the problem and do a single solve. havn't finished yet
function remake_single_solve(; model=model, signal=signal::Signal, prc2="NA", mute_parameter_disp=false)
    # ======== change parameter set with a particular prc2 rate
    p = vcat([collect(parameter_set[signal.db_idx, :]), signal.freq, 0.0, signal.phase]...)
    if prc2 != "NA"
        p[5] = prc2
    end
    pmap = parameters(model) .=> p
    u0 = collect(initial_condition[signal.db_idx, :])
    u0map = species(model) .=> u0
    reset_signal(; signal=signal, amplitude=signal.amplitude, freq=signal.freq)
    @show signal
    #  ======= construct prob and solve
    prob = ODEProblem(model, u0map, tspan, pmap)
    # prob = remake_prob(model, u0map, signal.tspan, p; prc2=prc2, mute_parameter_disp=false)
    sol = solve(prob, Rosenbrock23(), callback=cb, tstops=ts)
    return sol, ts
end



function single_solve_plot(; model=model, signal=signal::Signal, type="pulsatile", title="on")
    # ======= make adjustment for the signal
    # if signal.freq != 0
    #     signal.phase = 3*pi/2 - signal.freq*signal.T_init
    #     period = 2pi/signal.freq
    #     num_cycle= signal.ΔT_Dll4_ref/(period)
    #     if floor(num_cycle) < 1.0
    #         @show signal.ΔT = num_cycle* period
    #     else
    #         @show signal.ΔT =  floor(num_cycle)* period  - 0.01
    #     end 
    #     # signal.ΔT = floor(num_cycle)* period  - 0.01
    #     signal.tspan = (0, signal.T_init + signal.ΔT + signal.T_init)
    #     @show signal
    # else signal.freq == 0.0
    #     signal.phase = 0.0
    # end
    # dump(signal)
    # reset_signal(; signal = signal, amplitude=signal.amplitude, freq=signal.freq)
    # ======= make adjustment for the signal
    u0map, pmap, p, tspan, ts, cb, sol = single_solve(; model=model, signal=signal)
    # println("ts:",ts)
    @show pmap
    @show t_switching = switching_time(; sol=sol, pulse_period=signal.T_init:0.1:signal.T_init+signal.ΔT, idx=[6, 9], return_plot=false)
    t_switching = round(t_switching - signal.T_init, digits=2)
    if title == "on"
        plt = plot(sol,
            vars=[4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            title="A=$(signal.amplitude), freq=$(signal.freq), ST=$t_switching",
            titlefont=font(10, "Arial"),
            dpi=500)
    elseif title == "off"
        plt = plot(sol,
            vars=[4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            # title = "A : $amplitude, freq = $freq, switch time : $t_switching",
            dpi=500)
    end
    tt = ts[1]:0.01:ts[2]
    if signal.freq == 0
        plot!(plt, [0, ts[1], ts[2], tspan[end]], [0, 0, signal.amplitude, 0],
            label="Sustainable Input", seriestype=:steppre, line=(:dashdot, 2), alpha=0.8,
            # ylims = [0, 400],
            fill=(0, 0.3, :blue), color="black", dpi=300)
        return plt
    elseif signal.freq != 0
        if type == "pulsatile"
            pulse_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)))
            plot!(plt, tt, pulse_signal.(tt, signal.amplitude, signal.freq, signal.phase),
                label="Pulsatile Input", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
                # ylims = [0, 700],
                fill=(0, 0.3, :darkgreen), color="black", dpi=300)
            return plt
        elseif type == "bump"
            bump_signal(t, A, w, ϕ) = A * (abs(cos(w * t + ϕ)))
            plot!(plt, tt, bump_signal.(tt, signal.amplitude, signal.freq, signal.phase),
                label="Pulsatile Input", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
                # ylims = [0, 700],
                fill=(0, 0.3, :darkgreen), color="black", dpi=300)
            return plt
        end
    end
end


function remake_prob(model, u0map, tspan, p; prc2=0.4, mute_parameter_disp=false)
    if mute_parameter_disp == false
        @show p
        p[5] = prc2
        @show p
    elseif mute_parameter_disp == true
        p[5] = prc2
    end
    pmap = parameters(model) .=> p
    prob_new = ODEProblem(model, u0map, tspan, pmap)
end




# make animation
function anim_prc2_changing(range; model=model, u0=u0, tspan=tspan, p=p, cb=cb, ts=ts)
    anim = @animate for prc2 ∈ range
        prob = remake_prob(model, u0, tspan, p; prc2=prc2)
        @time sol = solve(prob, Rosenbrock23(), callback=cb, tstops=ts)
        plt = plot(sol,
            vars=[4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
            lw=1.5,
            xlabel="Time", ylabel="Concentration",
            title="PRC2 rate : $prc2",
            dpi=500)
    end
    return gif(anim, fps=1)
end


# ------ Find the switching time
function switching_time(; sol=sol, pulse_period=1.0:0.1:300, idx=[6, 9], return_plot=true)
    ∇_1(t, idx) = sol(t, Val{1}, idxs=idx)
    # plot(∇_1(pulse_period, idx))
    max_value, max_pos = findmax(∇_1(pulse_period, idx))
    t_max_id = max_pos[2]
    t_switching = pulse_period[t_max_id]
    # @show t_switching
    # sol(t_switching)
    # if return_plot == true
    #     plt_final = scatter!(plt, ones(2)*t_switching, sol(t_switching)[[6,9]], color = "purple", label = "Switching Event")
    #     display(plt_final)
    # end
    return t_switching
end





## =========================== reset signal =================
function reset_signal(; signal::Signal, tspan=(0.0, 600.0), amplitude=100, freq=0.1, T_init=100, ΔT_Dll4_ref=150)
    # signal.db_idx = db_idx # ---- no need for to state db_idx
    # signal.tspan = tspan
    signal.amplitude = amplitude
    signal.freq = freq
    signal.T_init = T_init
    signal.phase = 3 * pi / 2 - signal.freq * signal.T_init # 🔴reset phase
    signal.ΔT_Dll4_ref = ΔT_Dll4_ref
    if signal.freq == 0.0
        period = 2pi / 10^-5
    elseif signal.freq != 0.0
        period = 2pi / signal.freq
    end
    num_cycle = signal.ΔT_Dll4_ref / (period)
    if floor(num_cycle) < 1.0
        signal.ΔT = num_cycle * period
        println("ΔT is :\n", signal.ΔT)
    else
        signal.ΔT = floor(num_cycle) * period - 0.01
        println("ΔT is :\n", signal.ΔT)
    end
    signal.tspan = (0, signal.T_init + signal.ΔT + signal.T_init)
    return signal
end













"""
Input:
1. Amplitude range
2. Frequency range
3. Database number db_idx
Note: T_init, ΔT, tspan will take from the Environment, and they are printed out on screen.
# example
```jldoctest
A_ω_ϕ_st_relation(;amplitude_range = 100:20:300, freq_range = 0:0.01:0.4, db_idx = 301)
```
"""
function A_ω_ϕ_st_relation(; model=model, amplitude_range=100:20:300, freq_range=0:0.01:0.4, db_idx=301, T_init=100, ΔT_Dll4_ref=100, tspan_rt=4, prc2="NA", mute_parameter_disp=true)
    @show db_idx, prc2
    switch_amplitude = []
    switch_frequency = []
    switch_phase = []
    switch_time = []
    println("One can change the following time variables in the environment")
    tspan = (0.0, tspan_rt * ΔT_Dll4_ref)
    @show T_init, ΔT_Dll4_ref, tspan
    @showprogress for amplitude ∈ amplitude_range
        for freq_i ∈ freq_range # for each frequency, one need a re-phase the signal
            signal.db_idx = db_idx
            signal.tspan = tspan
            signal.amplitude = amplitude
            signal.freq = freq_i
            signal.phase = 3 * pi / 2 - signal.freq * signal.T_init
            signal.T_init = T_init
            signal.ΔT_Dll4_ref = ΔT_Dll4_ref
            period = 2pi / signal.freq
            num_cycle = signal.ΔT_Dll4_ref / (period)
            if floor(num_cycle) < 1.0
                @show signal.ΔT = num_cycle * period
            else
                @show signal.ΔT = floor(num_cycle) * period - 0.01
            end

            _, _, _, _, ts, _, single_sol = single_solve(; model=model, signal=signal::Signal, prc2=prc2, mute_parameter_disp=mute_parameter_disp)
            t_switching = switching_time(; sol=single_sol, pulse_period=signal.T_init:0.1:signal.T_init+signal.ΔT, idx=[6, 9], return_plot=false)
            check = check_switching(single_sol, ts, tspan)

            if check == -1
                push!(switch_time, t_switching)
                append!(switch_amplitude, amplitude)
                append!(switch_frequency, freq_i)
                # append!(switch_phase, phase)

                # plt = plot(single_sol,
                # # vars = [4, 5, 6, 9], # [MR,KDM5A,H4,H27]
                # vars = [ 6, 9], # [H4,H27]
                # lw = 2,
                # xlabel = "Time", ylabel = "Concentration",
                # foreground_color_legend = nothing,
                # # ylims = [0,500],
                # title = "A : $amplitude, freq = $freq_i, switch time : $t_switching",
                # dpi = 300)

            end

        end
    end
    A_ω_ϕ_st_database = DataFrame(freq=switch_frequency, amp=switch_amplitude, stime=switch_time)
    return A_ω_ϕ_st_database
end


# end

function df4d_sub_gen(df4d; fix_phase=true, fix_amp=true, amplitude_select=rand(unique(df4d.amp)))
    if fix_phase == true
        df4d_phase_0 = filter(row -> row.phase in [0], df4d)
        df4d_phase_1 = filter(row -> row.phase in [1], df4d)
        df4d_phase_2 = filter(row -> row.phase in [2], df4d)
        df4d_phase_3 = filter(row -> row.phase in [3], df4d)
        df4d_phase_4 = filter(row -> row.phase in [4], df4d)
        df4d_phase_5 = filter(row -> row.phase in [5], df4d)
        df4d_phase_6 = filter(row -> row.phase in [6], df4d)
    end
    if fix_amp == true
        df4d_amp_select = filter(row -> row.amp in [amplitude_select], df4d)
    end
    return df4d_phase_0, df4d_phase_1, df4d_phase_2, df4d_phase_3, df4d_phase_4, df4d_phase_5, df4d_phase_6, df4d_amp_select
end



function df4d_fix_phase_freq_vs_ST(df4d_phase_0; ΔT=100, amplitude_select=false, palette=:RdYlBu_6, save=false)
    @show db_idx
    if isempty(amplitude_select) == false
        df4d_phase_0_amp_select = filter(row -> row.amp in amplitude_select, df4d_phase_0)
        color_catg = length(amplitude_select)
    end
    plt = @df df4d_phase_0_amp_select plot(
        :freq,
        :stime,
        group=:amp,
        palette=palette,
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
    if save
        class = "positive_pulse" * "_id_$db_idx" * "/ΔT = $ΔT" * "/activation_time/"
        save_path = joinpath(@__DIR__, "figures", "switching_dynamics/$class") # generate path to save
        isdir(save_path) || mkpath(save_path) # create directory if not exist
        savefig(plt, save_path * "Fix_phase_|_ω_vs_ST" * ".png")
    end
    return plt
end






## ================== random process for signal =============================
# function poisson_distributed_puslses(freq, T_init; relax_mean = 10.0, ΔT_Dll4_ref = 150, num = 10, relax_to_period = true)
#     pulse_width = period = 2pi/freq # set the pulse width equals to the period
#     @show period
#     num = floor(Int, ΔT_Dll4_ref/period)
#     @show num

#     # set poisson intervals between pulses
#     if relax_to_period == true
#         relax_mean = pulse_width
#         poisson_interval = rand(Poisson(relax_mean), num)
#     elseif relax_to_period == false
#         poisson_interval = rand(Poisson(relax_mean), num)
#     end 


#     pulse_seq = []
#     for i in 1:length(poisson_interval)
#         push!(pulse_seq, pulse_width) # add pulse_width to pulse_seq
#         push!(pulse_seq, poisson_interval[i]) # add the current element to pulse_seq
#     end
#     @show pulse_seq
#     @show pulses_with_init = cumsum(pulse_seq)

#     return pulse_width, poisson_interval, pulse_seq, pulses_with_init
# end



function generate_poisson_pulses(freq, initial_time, max_time)
    pulses = [] # initialize the pulses with the initial time
    running_sum = initial_time # update the running sum
    pulse_width = period = 2pi / freq # set the pulse width equals to the period
    # @show period
    poisson_mean = pulse_width
    poisson = Poisson(poisson_mean) # create a Poisson distribution with the given mean

    # generate pulses and intervals until the running sum reaches max_time
    while running_sum < max_time
        if max_time - initial_time < period
            push!(pulses, [initial_time, max_time])
            break
        elseif max_time - initial_time >= period
            # generate a random interval from the Poisson distribution
            interval = rand(poisson)
            push!(pulses, running_sum .+ [0, pulse_width])
            # add the interval and the pulse_width to the running sum and append them to the pulses
            running_sum += pulse_width + interval
            # @show pulses
            # @show interval
            # @show running_sum
        end
    end
    if pulses[end][2] > max_time
        pulses = pulses[1:end-1]
    end
    return pulses
end


function heaviside(t)
    0.5 * (sign(t) + 1)
end

function interval(t, a, b)
    heaviside(t - a) - heaviside(t - b)
end

function compose_pulses(t, intervals)
    result = 0.0
    for i in 1:length(intervals)
        result += interval(t, intervals[i][1], intervals[i][2])
    end
    return result
end

"""
Model Description:
For a given frequency, `generate_poisson_pulses` will generate a instance of poisson pulses sequence
in the given duration of signal. Each model loading is associated with a fixed pulse sequence.
"""
function Load_poisson_model(freq, T_init, ΔT)
    single_instance_sequence = generate_poisson_pulses(freq, T_init, T_init + ΔT)
    # @show single_instance_sequence
    FR(t) = compose_pulses(t, single_instance_sequence) # each freq has a unique pattern
    # @register_symbolic FR(t)
    model_poisson_pulse = @reaction_network begin
        (A * FR(t), 1.0), R ↔ NR               # NICD binds RBPJ
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
    end k0 k1 k2 d m p k pp kk δ α1 A
    return model_poisson_pulse, single_instance_sequence
end


##  ======================================================= single solve for poisson pulses =============================
function single_solve_poisson(; signal=signal::Signal, prc2=prc2, mute_parameter_disp=mute_parameter_disp)
    db_idx = signal.db_idx
    tspan = signal.tspan
    T_init = signal.T_init
    ΔT = signal.ΔT
    freq = signal.freq
    amplitude = signal.amplitude

    # ====== each time loading poisson model will generate a sequence of poisson pulses within [T_init, T_init + ΔT]
    model, single_instance_sequence = Load_poisson_model(freq, T_init, ΔT)
    p = vcat([collect(parameter_set[db_idx, :]), 0.0]...)
    pmap = parameters(model) .=> p
    u0 = collect(initial_condition[db_idx, :])
    u0map = species(model) .=> u0
    ts, cb = make_cb([T_init, T_init + ΔT], 12, amplitude)
    prob1 = ODEProblem(model, u0map, tspan, pmap)
    if prc2 != "NA"
        prob1 = remake_prob(model, u0map, signal.tspan, p; prc2=prc2, mute_parameter_disp=mute_parameter_disp)
    end
    sol = solve(prob1, Rosenbrock23(), callback=cb, tstops=ts)
    return u0map, pmap, p, signal.tspan, ts, cb, sol
end





# 🍕 ======== generate the data set for poisson pulses ======
function A_ω_ϕ_st_relation_poisson(; amplitude_range=100:20:300, freq_range=0:0.01:0.4, db_idx=301, T_init=100, ΔT_Dll4_ref=100, tspan_rt=4, prc2="NA", mute_parameter_disp=true)
    @show db_idx, prc2
    switch_amplitude = []
    switch_frequency = []
    switch_phase = []
    switch_time = []
    println("One can change the following time variables in the environment")
    tspan = (0.0, tspan_rt * ΔT_Dll4_ref)
    @show T_init, ΔT_Dll4_ref, tspan
    @showprogress for amplitude ∈ amplitude_range
        for freq_i ∈ freq_range # for each frequency, one need a re-phase the signal
            signal.db_idx = db_idx
            signal.tspan = tspan
            signal.amplitude = amplitude
            signal.freq = freq_i
            signal.T_init = T_init

            # signal.ΔT_Dll4_ref = ΔT_Dll4_ref
            # period = 2pi/signal.freq
            # num_cycle= signal.ΔT_Dll4_ref/(period)
            # if floor(num_cycle) < 1.0
            #     @show signal.ΔT = num_cycle* period
            # else
            #     @show signal.ΔT =  floor(num_cycle)* period  - 0.01
            # end 

            _, _, _, _, ts, _, single_sol_poisson = single_solve_poisson(; signal=signal::Signal, prc2=prc2, mute_parameter_disp=mute_parameter_disp)
            t_switching = switching_time(; sol=single_sol_poisson, pulse_period=signal.T_init:0.1:signal.T_init+signal.ΔT, idx=[6, 9], return_plot=false)
            check = check_switching(single_sol_poisson, ts, tspan)

            if check == -1
                push!(switch_time, t_switching)
                append!(switch_amplitude, amplitude)
                append!(switch_frequency, freq_i)
                # append!(switch_phase, phase)

                # plt = plot(single_sol,
                # # vars = [4, 5, 6, 9], # [MR,KDM5A,H4,H27]
                # vars = [ 6, 9], # [H4,H27]
                # lw = 2,
                # xlabel = "Time", ylabel = "Concentration",
                # foreground_color_legend = nothing,
                # # ylims = [0,500],
                # title = "A : $amplitude, freq = $freq_i, switch time : $t_switching",
                # dpi = 300)

            end

        end
    end
    A_ω_ϕ_st_database_poisson = DataFrame(freq=switch_frequency, amp=switch_amplitude, stime=switch_time)
    return A_ω_ϕ_st_database_poisson
end



# ======== two genes comparison and prc comparison =================================
function Two_Genes_TS_by_Prc2(; model=model, signal1, signal2, prc2=0.1, tspan_rt=4, title_on=true, legend_title_on=true, vars_to_show=[4, 5, 6, 9], type="pulsatile")
    # ======== Gene 1 with Dll4 sustainable signal
    # sol_gene1, ts1 = remake_prob(; model = model, amplitude = amplitude1, T_init = T_init, ΔT = ΔT, tspan_rt = tspan_rt, prc2 = prc2)
    p = vcat([collect(parameter_set[signal1.db_idx, :]), signal1.freq, 0.0, signal1.phase]...)
    # @show p
    u0 = collect(initial_condition[signal1.db_idx, :])
    u0map = species(model) .=> u0
    @show u0map
    # ======= since the gene1 part is for Dll4, so we set the signal1.freq to 0.0
    reset_signal(; signal=signal1, amplitude=signal1.amplitude, freq=0)
    @show prob1 = remake_prob(model, u0map, signal1.tspan, p; prc2=prc2, mute_parameter_disp=false)
    u0map, pmap, p, tspan, ts1, cb, sol_gene1 = single_solve(; model=model_pulsatile, signal=signal1)

    @show t_switching1 = switching_time(; sol=sol_gene1, pulse_period=signal1.T_init:0.1:signal1.T_init+signal1.ΔT, idx=[6, 9], return_plot=false)
    t_switching1 = t_switching1 - signal1.T_init

    if title_on == true && legend_title_on == true
        plt_gene1_Dll4 = plot(sol_gene1,
            vars=vars_to_show, # [MR,KDM5A,H4,H27,KDM6A]
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            legend_title="Gene 1", legendtitlefontsize=10,
            title="A=$(signal1.amplitude), freq=0, PRC2 rate=$prc2,  ST=$t_switching1",
            titlefont=font(10, "Arial"),
            dpi=500)
    elseif title_on == false && legend_title_on == true
        plt_gene1_Dll4 = plot(sol_gene1,
            vars=vars_to_show, # [MR,KDM5A,H4,H27,KDM6A]
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            legend_title="Gene 1", legendtitlefontsize=10,
            # title = "A=$amplitude1, freq=0, PRC2 rate=$prc2,  ST=$t_switching1",
            titlefont=font(10, "Arial"),
            dpi=500)
    elseif title_on == true && legend_title_on == false
        plt_gene1_Dll4 = plot(sol_gene1,
            vars=vars_to_show, # [MR,KDM5A,H4,H27,KDM6A]
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            # legend_title = "Gene 1", legendtitlefontsize = 10,
            title="A=$(signal1.amplitude), freq=0, PRC2 rate=$prc2,  ST=$t_switching1",
            titlefont=font(10, "Arial"),
            dpi=500)
    elseif title_on == false && legend_title_on == false
        plt_gene1_Dll4 = plot(sol_gene1,
            vars=vars_to_show, # [MR,KDM5A,H4,H27,KDM6A]
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            # legend_title = "Gene 1", legendtitlefontsize = 10,
            # title = "A=$amplitude1, freq=0, PRC2 rate=$prc2,  ST=$t_switching1",
            titlefont=font(10, "Arial"),
            dpi=500)
    end
    plot!(plt_gene1_Dll4, [0, ts1[1], ts1[2], tspan_rt * signal1.ΔT[end]], [0, 0, signal1.amplitude, 0],
        label="Sustained Input", seriestype=:steppre, line=(:dashdot, 2), alpha=0.8,
        # ylims = [0, 400],
        fill=(0, 0.3, :blue), color="black", dpi=500)



    # ======= Gene 2 with Dll1 pulsatile signal
    osci_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)))
    p = vcat([collect(parameter_set[signal2.db_idx, :]), signal2.freq, 0.0, signal2.phase]...)
    # @show p
    u0 = collect(initial_condition[signal2.db_idx, :])
    u0map = species(model) .=> u0
    @show u0map
    # ======= since the gene1 part is for Dll1, so we let the frequency to be the signal2.freq 
    reset_signal(; signal=signal2, amplitude=signal2.amplitude, freq=signal2.freq)
    @show prob2 = remake_prob(model, u0map, signal2.tspan, p; prc2=prc2, mute_parameter_disp=false)
    u0map, pmap, p, tspan, ts1, cb, sol_gene2 = single_solve(; model=model_pulsatile, signal=signal2)

    @show t_switching2 = switching_time(; sol=sol_gene2, pulse_period=signal2.T_init:0.1:signal2.T_init+signal2.ΔT, idx=[6, 9], return_plot=false)
    t_switching2 = t_switching2 - signal2.T_init

    # sol_gene2, ts2 = remake_solve_prob_by_ID(; model=model, db_idx=id2, freq=id2_freq, phase=phase2, amplitude=amplitude2, T_init=T_init, ΔT=ΔT, tspan_rt=tspan_rt, prc2=prc2)
    # @show t_switching2 = switching_time(; sol=sol_gene2, pulse_period=T_init:0.1:T_init+ΔT, idx=[6, 9], return_plot=false)
    # t_switching2 = t_switching2 - T_init


    tt = ts2[1]:0.01:ts2[2]
    if title_on == true && legend_title_on == true
        plt_gene2_Dll1 = plot(sol_gene2,
            vars=vars_to_show, # [MR,KDM5A,H4,H27,KDM6A]
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            legend_title="Gene 2", legendtitlefontsize=10,
            title="A=$amplitude2, freq=$id2_freq, PRC2 rate=$prc2, ST=$t_switching2",
            titlefont=font(10, "Arial"),
            dpi=500)
    elseif title_on == false && legend_title_on == true
        plt_gene2_Dll1 = plot(sol_gene2,
            vars=vars_to_show, # [MR,KDM5A,H4,H27,KDM6A]
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            legend_title="Gene 2", legendtitlefontsize=10,
            # title = "A=$amplitude2, freq=$id2_freq, PRC2 rate=$prc2, ST=$t_switching2",
            titlefont=font(10, "Arial"),
            dpi=500)
    elseif title_on == true && legend_title_on == false
        plt_gene2_Dll1 = plot(sol_gene2,
            vars=vars_to_show, # [MR,KDM5A,H4,H27,KDM6A]
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            # legend_title = "Gene 2", legendtitlefontsize = 10,
            title="A=$amplitude2, freq=$id2_freq, PRC2 rate=$prc2, ST=$t_switching2",
            titlefont=font(10, "Arial"),
            dpi=500)
    elseif title_on == false && legend_title_on == false
        plt_gene2_Dll1 = plot(sol_gene2,
            vars=vars_to_show, # [MR,KDM5A,H4,H27,KDM6A]
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            # legend_title = "Gene 1", legendtitlefontsize = 10,
            # title = "A=$amplitude2, freq=$id2_freq, PRC2 rate=$prc2, ST=$t_switching2",
            titlefont=font(10, "Arial"),
            dpi=500)
    end

    # plot!(plt_gene2_Dll1, tt, osci_signal.(tt, amplitude2, id2_freq, phase2),
    #     label = "Pulsatile Input", seriestype = :steppre, line = (:dot, 2), alpha = 0.8,
    #     # ylims = [0, 700],
    #     fill = (0, 0.3, :darkgreen), color = "black", dpi = 500)

    if id2_freq != 0
        if type == "pulsatile"
            pulse_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)))
            plot!(plt_gene2_Dll1, tt, pulse_signal.(tt, amplitude2, id2_freq, phase2),
                label="Pulsatile Input", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
                # ylims = [0, 700],
                fill=(0, 0.3, :darkgreen), color="black", dpi=300)
            # return plt_gene2_Dll1
        elseif type == "bump"
            bump_signal(t, A, w, ϕ) = A * (abs(cos(w * t + ϕ)))
            plot!(plt_gene2_Dll1, tt, bump_signal.(tt, amplitude2, id2_freq, phase2),
                label="Pulsatile Input", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
                # ylims = [0, 700],
                fill=(0, 0.3, :darkgreen), color="black", dpi=300)
            # return plt_gene2_Dll1
        end
    end


    plt_combine = plot(plt_gene1_Dll4, plt_gene2_Dll1, layout=(2, 1))
    return plt_gene1_Dll4, plt_gene2_Dll1, plt_combine
end
