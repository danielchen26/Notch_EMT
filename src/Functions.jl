## ====================== Defining the model and parameters ==========================
module model_parameters
Base.@kwdef mutable struct Signal
    db_idx::Int64
    tspan::Tuple{Float64,Float64}
    freq::Float64
    phase::Float64
    amplitude::Float64
    T_init::Float64
    ΔT::Float64
end
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
        end k0 k1 k2 d m p k pp kk δ α1 w A ϕ # put A at last as the control 13th variable
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
# add document for the following function
# ts_in is the callback function working time duration, ts_in[1] is the start time, ts_in[2] is the end time
# at the start time, the value for 13th paraemter is A, at the end time, the value for 13th parameter is turned off to 0.
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

#* an alternative way to make callbacks, same as above
"""
# example
```julia
ps_cb_1 = PresetTimeCallback([T_init], integ -> integ.p[end-1] = amplitude)
ps_cb_2 = PresetTimeCallback([T_init + ΔT], integ -> integ.p[end-1] = 0.0)

sol = solve(deepcopy(prob1); callback=CallbackSet(ps_cb_1, ps_cb_2))
plot(sol, vars=[4, 5, 6, 9], lw=1.5, xlabel="Time", ylabel="Concentration", dpi=500)
```
"""
function make_callbacks(T_init, ΔT, amplitude)
    ps_cb_1 = PresetTimeCallback([T_init], integ -> integ.p[end-1] = amplitude)
    ps_cb_2 = PresetTimeCallback([T_init + ΔT], integ -> integ.p[end-1] = 0.0)
    return callback = CallbackSet(ps_cb_1, ps_cb_2)
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




# write a function to extract for each frequency the minimum driving amplitude
function extract_min_amp(df::DataFrame)
    # Group the dataframe by the 'freq' column
    grouped_df = groupby(df, :freq)

    # Initialize an empty dataframe to store the results
    min_amp_df = DataFrame(freq=Float64[], min_amp=Float64[])

    # Iterate over the groups and find the minimum amp for each freq
    for group in grouped_df
        min_amp = minimum(group[:, :amp])
        freq = group[1, :freq]
        push!(min_amp_df, (freq, min_amp))
    end

    return rename!(min_amp_df, :min_amp => :amp)
end


function add_Boundary_event(switch_amplitude, switch_frequency, plt)
    df = DataFrame(switch_amplitude=switch_amplitude, switch_frequency=switch_frequency)
    gp = groupby(df, :switch_frequency)
    boundary_amplitude = []
    for i in gp
        append!(boundary_amplitude, minimum(i.switch_amplitude))
    end
    critical_freq = [keys(gp)[i].switch_frequency for i = 1:gp.ngroups]
    plot!(plt, boundary_amplitude, critical_freq, 
        label="Switching Boundary", 
        linewidth=3,
        color=:red)
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
function single_solve(; model=model, db_idx, freq, phase, amplitude, T_init, ΔT, tspan, prc2="NA", mute_parameter_disp=false, phase_reset=true) #🍏🔴added prc2 changing option
    if phase_reset == true && freq != 0.0
        phase = 3pi / 2 - freq * T_init
    end
    p = vcat([collect(parameter_set[db_idx, :]), freq, 0.0, phase]...)
    pmap = parameters(model) .=> p
    u0 = collect(initial_condition[db_idx, :])
    u0map = species(model) .=> u0
    ts, cb = make_cb([T_init, T_init + ΔT], 13, amplitude)
    prob1 = ODEProblem(model, u0map, tspan, pmap)
    if prc2 != "NA"
        prob1 = remake_prob(model, u0map, tspan, p; prc2=prc2, mute_parameter_disp=mute_parameter_disp)
    end
    sol = solve(prob1, Rosenbrock23(), callback=cb, tstops=ts)
    return u0map, pmap, p, tspan, ts, cb, sol, phase
end


function single_solve_plot(; model=model, db_idx, phase, freq, amplitude, T_init, ΔT, type="pulsatile", title="on", phase_reset=true, T_final = 3*ΔT, arg...)
    tspan = (0.0, T_init + T_final)
    u0map, pmap, p, tspan, ts, cb, sol, phase = single_solve(; model=model, db_idx=db_idx, freq=freq, phase=phase, amplitude=amplitude, T_init=T_init, ΔT=ΔT, tspan=tspan, phase_reset=phase_reset, arg...)
    # @show pmap
    check = check_switching(sol, ts, tspan)
    println(check == -1 ? "Histone states swithced" : "NO swithing")
    # t_switching = switching_time(; sol=sol, pulse_period=T_init:0.1:T_init+ΔT, idx=[6, 9], return_plot=false)
    # t_switching = round(t_switching - T_init, digits=2)
    t_switching = isnothing(find_ST(sol)) ? "No Swithing" : round(find_ST(sol), digits=2)
    @show t_switching
    if title == "on"
        plt = plot(sol,
            vars=[4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
            # vars = [1,2,3,4,5,6,7,8,9,10,11], # all variables
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            dpi=500)
        check == -1 ?
        title!(plt, " A=$amplitude, freq=$freq, ST=$t_switching",
            titlefont=font(10, "Arial")) :
        title!(plt, " A=$amplitude, freq=$freq, ST= No Swithing",
            titlefont=font(10, "Arial"))
    elseif title == "off"
        plt = plot(sol,
            vars=[4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            dpi=500)
    end
    tt = ts[1]:0.01:ts[2]
    if freq == 0
        plot!(plt, [0, ts[1], ts[2], tspan[end]], [0, 0, 2 * amplitude, 0], # FIXME: need to change to 2*amplitude
            label="Sustained Input", seriestype=:steppre, line=(:dashdot, 2), alpha=0.8,
            # ylims = [0, 400],
            fill=(0, 0.3, :blue), color="black", dpi=300)
        return plt
    elseif freq != 0
        if type == "pulsatile"
            pulse_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)))
            plot!(plt, tt, pulse_signal.(tt, amplitude, freq, phase),
                label="Pulsatile Input", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
                # ylims = [0, 700],
                fill=(0, 0.3, :darkgreen), color="black", dpi=300)
            return plt
        elseif type == "bump"
            bump_signal(t, A, w, ϕ) = A * (abs(cos(w * t + ϕ)))
            plot!(plt, tt, bump_signal.(tt, amplitude, freq, phase),
                label="Bump Input", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
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

"""
Input:
1. Database ID : `db_idx`
2. Frequency
3. Phase
4. Amplitude
5. T_init
6. ΔT
7. tspan
8. prc2 rate : `specify a new one`
# example
```jldoctest
sol = remake_solve_prob_by_ID(model, db_idx, freq, phase, amplitude, T_init, ΔT, tspan ; prc2 = 0.4)
```
"""
function remake_solve_prob_by_ID(; model, db_idx, freq=0.0, phase=0.0, amplitude, T_init=100, ΔT=100, tspan_rt=4, prc2=0.4)
    tspan = (0.0, tspan_rt * ΔT)
    @show freq
    @show phase
    @show T_init, ΔT, tspan
    println("\n")
    # ======== change parameter set with a particular prc2 rate
    p = vcat([collect(parameter_set[db_idx, :]), freq, 0.0, phase]...)
    println("Parameter Set ID $db_idx : ", p)
    p[5] = prc2
    println("New PRC2 rate for ID $db_idx : ", p)
    pmap = parameters(model) .=> p
    u0 = collect(initial_condition[db_idx, :])
    u0map = species(model) .=> u0
    ts, cb = make_cb([T_init, T_init + ΔT], 13, amplitude)
    #  ======= construct prob and solve
    prob = ODEProblem(model, u0map, tspan, pmap)
    sol = solve(prob, Rosenbrock23(), callback=cb, tstops=ts)
    return sol, ts, tspan
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
#! this function might find switching time better 
function find_ST(sol::ODESolution)
    t = sol.t
    y1 = sol[6, :] #H4
    y2 = sol[9, :] #27
    for i in 1:length(t)-1
        if (y1[i] - y2[i]) * (y1[i+1] - y2[i+1]) < 0
            # trajectories cross between i and i+1
            return (t[i] + t[i+1]) / 2
        end
    end
    return nothing  # trajectories never cross
end


#! fixed the phase2 for the pulsalatile signle for gene 2
"""
    Two_Genes_TS_by_Prc2(; model = model, id1 = 592, id2 = 49, id2_freq = 0.2, phase2 = 0.0, amplitude1 = 130, amplitude2 = 130, prc2 = 0.1, T_init = 100, ΔT = 100, tspan_rt = 4, title_on = true, legend_title_on = true, vars_to_show = [4, 5, 6, 9], type = "pulsatile", phase2_reset = true )

TBW
"""
function Two_Genes_TS_by_Prc2(; model=model, id1=592, id2=49, id2_freq=0.2, phase2=0.0, amplitude1=130, amplitude2=130, prc2=0.1, T_init=1e-10, ΔT=100, tspan_rt=2, title_on=true, legend_title_on=true, vars_to_show=[4, 5, 6, 9], type="pulsatile", phase2_reset=true)
    # ======== signal form =================
    osci_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)))
    # ======== Gene 1 with Dll4 sustainable signal
    sol_gene1, ts1, tspan1 = remake_solve_prob_by_ID(; model=model, db_idx=id1, amplitude=amplitude1, T_init=T_init, ΔT=ΔT, tspan_rt=tspan_rt, prc2=prc2) #! default freq = 0.0, phase = 0.0 which associates with Dll4
    check1 = check_switching(sol_gene1, ts1, tspan1)
    t_switching1 = switching_time(; sol=sol_gene1, pulse_period=T_init:0.1:T_init+ΔT, idx=[6, 9], return_plot=false)
    t_switching1 = round(t_switching1 - T_init, digits=3)
    t_switching1_improved = isnothing(find_ST(sol_gene1)) ? "No Switching" : round(find_ST(sol_gene1), digits=2)
    @show t_switching1, t_switching1_improved


    if title_on == true && legend_title_on == true
        plt_gene1_Dll4 = plot(sol_gene1,
            vars=vars_to_show, # [MR,KDM5A,H4,H27,KDM6A]
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            legend_title="Gene 1", legendtitlefontsize=10,
            dpi=500)
        check1 == -1 ?
        title!(plt_gene1_Dll4, "A=$(osci_signal(0,amplitude1,0,0)), freq=0, PRC2 rate=$prc2,  ST=$t_switching1_improved",
            titlefont=font(10, "Arial")) :
        title!(plt_gene1_Dll4, "A=$(osci_signal(0,amplitude1,0,0)), freq=0, PRC2 rate=$prc2,  ST=No Switching",
            titlefont=font(10, "Arial"))
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
            dpi=500)
        check1 == -1 ?
        title!(plt_gene1_Dll4, "A=$(osci_signal(0,amplitude1,0,0)), freq=0, PRC2 rate=$prc2,  ST=$t_switching1_improved",
            titlefont=font(10, "Arial"),) :
        title!(plt_gene1_Dll4, "A=$(osci_signal(0,amplitude1,0,0)), freq=0, PRC2 rate=$prc2,  ST=No Switching",
            titlefont=font(10, "Arial"))
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
    plot!(plt_gene1_Dll4, [0, ts1[1], ts1[2], tspan_rt * ΔT[end]], [0, 0, osci_signal(0, amplitude1, 0, 0), 0],
        label="Sustained Signal", seriestype=:steppre, line=(:dashdot, 2), alpha=0.8,
        # ylims = [0, 400],
        fill=(0, 0.3, :blue), color="black", dpi=500)

    # * =======================  Gene 2 with Dll1 pulsatile signal =============================
    # osci_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)))
    if phase2_reset == true && id2_freq != 0 #! phase should not be reseted for Dll4 sustained signal
        phase2 = 3pi / 2 - id2_freq * T_init #! fixed the phase2 for the pulsatile signle for gene 2
    end
    sol_gene2, ts2, tspan2 = remake_solve_prob_by_ID(; model=model, db_idx=id2, freq=id2_freq, phase=phase2, amplitude=amplitude2, T_init=T_init, ΔT=ΔT, tspan_rt=tspan_rt, prc2=prc2)
    check2 = check_switching(sol_gene2, ts2, tspan2)
    # t_switching2 = switching_time(; sol=sol_gene2, pulse_period=T_init:0.1:T_init+ΔT, idx=[6, 9], return_plot=false)
    # t_switching2 = round(t_switching2 - T_init, digits=3)
    t_switching2_improved = isnothing(find_ST(sol_gene2)) ? "No Switching" : round(find_ST(sol_gene2), digits=2)
    @show t_switching2_improved

    tt = ts2[1]:0.01:ts2[2]
    if title_on == true && legend_title_on == true
        plt_gene2_Dll1 = plot(sol_gene2,
            vars=vars_to_show, # [MR,KDM5A,H4,H27,KDM6A]
            lw=2,
            xlabel="Time", ylabel="Concentration",
            foreground_color_legend=nothing,
            legend_title="Gene 2", legendtitlefontsize=10,
            dpi=500)
        check2 == -1 ?
        title!(plt_gene2_Dll1, "A=$(osci_signal(0,amplitude2,id2_freq,3pi/2 + 1e-10)), freq=$id2_freq, PRC2 rate=$prc2,  ST=$t_switching2_improved",
            titlefont=font(10, "Arial"),) :
        title!(plt_gene2_Dll1, "A=$(osci_signal(0,amplitude2,id2_freq,3pi/2 + 1e-10)), freq=$id2_freq, PRC2 rate=$prc2,  ST=No Switching",
            titlefont=font(10, "Arial"))
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
            dpi=500)
        check2 == -1 ?
        title!(plt_gene2_Dll1, "A=$(osci_signal(0,amplitude2,id2_freq,3pi/2 + 1e-10 )), freq=$id2_freq, PRC2 rate=$prc2,  ST=$t_switching2_improved",
            titlefont=font(10, "Arial"),) :
        title!(plt_gene2_Dll1, "A=$(osci_signal(0,amplitude2,id2_freq,3pi/2 + 1e-10 )), freq=$id2_freq, PRC2 rate=$prc2,  ST=No Switching",
            titlefont=font(10, "Arial"))
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
                label="Pulsatile Signal", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
                # ylims = [0, 700],
                fill=(0, 0.3, :darkgreen), color="black", dpi=300)
            # return plt_gene2_Dll1
        elseif type == "bump"
            bump_signal(t, A, w, ϕ) = A * (abs(cos(w * t + ϕ)))
            plot!(plt_gene2_Dll1, tt, bump_signal.(tt, amplitude2, id2_freq, phase2),
                label="Pulsatile Signal", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
                # ylims = [0, 700],
                fill=(0, 0.3, :darkgreen), color="black", dpi=300)
            check2 == -1 ?
            title!(plt_gene2_Dll1, "A=$amplitude2, freq=$id2_freq, PRC2 rate=$prc2,  ST=$t_switching2_improved",
                titlefont=font(10, "Arial"),) :
            title!(plt_gene2_Dll1, "A=$amplitude2, freq=$id2_freq, PRC2 rate=$prc2,  ST=No Switching",
                titlefont=font(10, "Arial"))
            # return plt_gene2_Dll1
        end
    end


    plt_combine = plot(plt_gene1_Dll4, plt_gene2_Dll1, layout=(2, 1))
    return plt_gene1_Dll4, plt_gene2_Dll1, plt_combine
end



function find_id_Dll4_vs_Dll1(id; amplitude=65, freq=0.15, ΔT=100, prc2=0.41, type="pulsatile")
    @show id
    single_gene_id = id
    id2_freq = freq
    phase2 = 0
    amplitude1 = amplitude2 = amplitude
    T_init = 1e-10
    ΔT = ΔT
    prc2 = prc2

    plt_gene2_Dll4, plt_gene2_Dll1, plt_2genes_compare_id =
        Two_Genes_TS_by_Prc2(; model=model_pulsatile,
            id1=single_gene_id, id2=single_gene_id,
            id2_freq=id2_freq, phase2=phase2, amplitude1=amplitude1,
            amplitude2=amplitude2, prc2=prc2,
            T_init=T_init, ΔT=ΔT, title_on=false, legend_title_on=false,
            type=type)
    title!(plt_2genes_compare_id, "id:$id")
    display(plt_2genes_compare_id)
end



# generate an animation for Two_Genes_TS_by_Prc2() for various prc2 values

function Two_genes_prc2_TS_animation(; prc2_range=0:0.02:1, model=model, id1, id2, id2_freq=0.2, phase2=0.0, amplitude1=130, amplitude2=130, T_init=1e-10, ΔT=100, tspan_rt=2, type="pulsatile", phase2_reset=true)
    osci_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)))
    anim = @animate for prc2 ∈ prc2_range
        # ======== Gene 1 with Dll4 sustainable signal
        _, _, plt_combine = Two_Genes_TS_by_Prc2(; model=model, id1=single_gene_id, id2=single_gene_id,
            id2_freq=id2_freq, phase2=phase2, amplitude1=amplitude1,
            amplitude2=amplitude2, prc2=prc2,
            T_init=T_init, ΔT=ΔT, title_on=true, legend_title_on=false,
            vars_to_show=[5, 6, 9], #tspan_rt = 2, # 
            type="pulsatile",
            phase2_reset=true)
    end
    return anim
end





# function Two_genes_prc2_TS_animation(; prc2_range=0:0.02:1, model=model, db_idx1, db_idx2, id2_freq=0.2, phase2=0.0, amplitude1=130, amplitude2=130, T_init=1e-10, ΔT=100, tspan_rt=2, type="pulsatile", phase2_reset=true)
#     osci_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)))
#     anim = @animate for prc2 ∈ prc2_range
#         # ======== Gene 1 with Dll4 sustainable signal
#         @time sol_gene1, ts1, tspan1 = remake_solve_prob_by_ID(; model=model, db_idx=db_idx1, amplitude=amplitude1, T_init=T_init, ΔT=ΔT, tspan_rt=tspan_rt, prc2=prc2)
#         check1 = check_switching(sol_gene1, ts1, tspan1)
#         t_switching1 = switching_time(; sol=sol_gene1, pulse_period=T_init:0.1:T_init+ΔT, idx=[6, 9], return_plot=false)
#         t_switching1 = round(t_switching1 - T_init, digits=2)
#         @show t_switching1

#         plt_gene1_Dll4 = plot(sol_gene1,
#             vars=[4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
#             lw=2,
#             xlabel="Time", ylabel="Concentration",
#             foreground_color_legend=nothing,
#             legend_title="Gene 1", legendtitlefontsize=10,
#             title="A=$(osci_signal(0,amplitude1,0,0)), freq=0, PRC2 rate=$prc2,  ST=$t_switching1",
#             titlefont=font(10, "Arial"),
#             dpi=500)
#         plot!(plt_gene1_Dll4, [0, ts1[1], ts1[2], tspan_rt * ΔT], [0, 0, osci_signal(0, amplitude1, 0, 0), 0],
#             label="Sustained Signal", seriestype=:steppre, line=(:dashdot, 2), alpha=0.8,
#             # ylims = [0, 400],
#             fill=(0, 0.3, :blue), color="black", dpi=500)

#         # ======= Gene 2 with Dll1 pulsatile signal
#         if phase2_reset == true && id2_freq != 0 #! phase should not be reseted for Dll4 sustained signal
#             phase2 = 3pi / 2 - id2_freq * T_init #! fixed the phase2 for the pulsatile signle for gene 2
#         end
#         sol_gene2, ts2, tspan2 = remake_solve_prob_by_ID(; model=model, db_idx=db_idx2, freq=id2_freq, phase=phase2, amplitude=amplitude2, T_init=T_init, ΔT=ΔT, tspan_rt=tspan_rt, prc2=prc2)
#         check2 = check_switching(sol_gene2, ts2, tspan2)
#         t_switching2 = switching_time(; sol=sol_gene2, pulse_period=T_init:0.1:T_init+ΔT, idx=[6, 9], return_plot=false)
#         t_switching2 = round(t_switching2 - T_init, digits=2)
#         @show t_switching2

#         tt = ts2[1]:0.01:ts2[2]
#         plt_gene2_Dll1 = plot(sol_gene2,
#             vars=[4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
#             lw=2,
#             xlabel="Time", ylabel="Concentration",
#             foreground_color_legend=nothing,
#             legend_title="Gene 2", legendtitlefontsize=10,
#             title="A=$(osci_signal(0,amplitude2,id2_freq,3pi/2 + 1e-10)), freq=$id2_freq, PRC2 rate=$prc2, ST=$t_switching2",
#             titlefont=font(10, "Arial"),
#             dpi=500)

#         if id2_freq != 0
#             if type == "pulsatile"
#                 pulse_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)))
#                 plot!(plt_gene2_Dll1, tt, pulse_signal.(tt, amplitude2, id2_freq, phase2),
#                     label="Pulsatile Signal", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
#                     # ylims = [0, 700],
#                     fill=(0, 0.3, :darkgreen), color="black", dpi=500)
#                 # return plt_gene2_Dll1
#             elseif type == "bump"
#                 bump_signal(t, A, w, ϕ) = A * (abs(cos(w * t + ϕ)))
#                 plot!(plt_gene2_Dll1, tt, bump_signal.(tt, amplitude2, id2_freq, phase2),
#                     label="Pulsatile Input", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
#                     # ylims = [0, 700],
#                     fill=(0, 0.3, :darkgreen), color="black", dpi=500)
#                 # return plt_gene2_Dll1
#             end
#         end

#         # plot!(plt_gene2_Dll1, tt, osci_signal.(tt, amplitude2, id2_freq, 0.0),
#         #     label = "Pulsatile Input", seriestype = :steppre, line = (:dot, 2), alpha = 0.8,
#         #     # ylims = [0, 700],
#         #     fill = (0, 0.3, :darkgreen), color = "black", dpi = 500)

#         # @time sol_gene2  = remake_solve_prob_by_ID(;model = model, db_idx = db_idx2, freq = 0.2,  amplitude = 130, T_init = T_init, ΔT = ΔT, tspan_rt = tspan_rt, prc2 = prc2)
#         # plt_gene2 = plot(sol_gene2,
#         #     vars = [4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
#         #     lw = 1.5,
#         #     xlabel = "Time", ylabel = "Concentration",
#         #     title = "PRC2 rate : $prc2",
#         #     dpi = 500)
#         plt_combine = plot(plt_gene1_Dll4, plt_gene2_Dll1, layout=(2, 1))
#     end
#     return gif(anim, fps=1)
# end

















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
function A_ω_st_relation(; model=model, db_idx=301, amplitude_range=100:20:300, freq_range=0:0.01:0.4, T_init=1e-10, ΔT=100, tspan_rt=2, prc2="NA", mute_parameter_disp=true)
    @show db_idx, prc2
    switch_amplitude = []
    switch_frequency = []
    switch_time = []
    switch_prc2 = []
    println("One can change the following time variables in the environment")
    tspan = (0.0, tspan_rt * ΔT)
    @show T_init, ΔT, tspan
    @showprogress for amplitude ∈ amplitude_range
        for freq_i ∈ freq_range # 
            freq_i != 0 ? phase = 3pi / 2 - freq_i * T_init : phase = 0 # reset phase when freq != 0
            u0map, pmap, p, tspan, ts, cb, single_sol, _ = single_solve(; model=model, db_idx=db_idx, freq=freq_i, phase=phase, amplitude=amplitude, T_init=T_init, ΔT=ΔT, tspan=tspan, prc2=prc2, mute_parameter_disp=mute_parameter_disp)
            # t_switching = switching_time(; sol=single_sol, pulse_period=T_init:0.1:T_init+ΔT, idx=[6, 9], return_plot=false)
            t_switching_improved = isnothing(find_ST(single_sol)) ? "No Switching" : round(find_ST(single_sol), digits=2)
            check = check_switching(single_sol, ts, tspan)
            if check == -1
                push!(switch_time, t_switching_improved)
                append!(switch_amplitude, amplitude)
                append!(switch_frequency, freq_i)
                append!(switch_prc2, p[5])
            end
        end
    end
    A_ω_st_database = DataFrame(freq=switch_frequency, amp=switch_amplitude, stime=switch_time, prc2=switch_prc2)
    return A_ω_st_database
end


function A_ω_st_relation_prc2_range(; model=model, db_idx=301, amplitude_range=100:20:300, freq_range=0:0.01:0.4, T_init=1e-10, ΔT=100, tspan_rt=2, prc2_range="NA", mute_parameter_disp=true)
    all_dataframes = []

    for prc2 in prc2_range
        @show db_idx, prc2
        switch_amplitude = []
        switch_frequency = []
        switch_time = []
        switch_prc2 = []
        println("One can change the following time variables in the environment")
        tspan = (0.0, tspan_rt * ΔT)
        @show T_init, ΔT, tspan
        @showprogress for amplitude ∈ amplitude_range
            for freq_i ∈ freq_range
                freq_i != 0 ? phase = 3pi / 2 - freq_i * T_init : phase = 0
                u0map, pmap, p, tspan, ts, cb, single_sol, _ = single_solve(; model=model, db_idx=db_idx, freq=freq_i, phase=phase, amplitude=amplitude, T_init=T_init, ΔT=ΔT, tspan=tspan, prc2=prc2, mute_parameter_disp=mute_parameter_disp)
                t_switching_improved = isnothing(find_ST(single_sol)) ? "No Switching" : round(find_ST(single_sol), digits=2)
                check = check_switching(single_sol, ts, tspan)
                if check == -1
                    push!(switch_time, t_switching_improved)
                    append!(switch_amplitude, amplitude)
                    append!(switch_frequency, freq_i)
                    append!(switch_prc2, p[5])
                end
            end
        end
        A_ω_st_database = DataFrame(freq=switch_frequency, amp=switch_amplitude, stime=switch_time, prc2=switch_prc2)
        push!(all_dataframes, A_ω_st_database)
    end

    combined_dataframe = vcat(all_dataframes...)
    return combined_dataframe
end






# function df4d_sub_gen(df4d; fix_phase=true, fix_amp=true, amplitude_select=rand(unique(df4d.amp)))
#     if fix_phase == true
#         df4d_phase_0 = filter(row -> row.phase in [0], df4d)
#         df4d_phase_1 = filter(row -> row.phase in [1], df4d)
#         df4d_phase_2 = filter(row -> row.phase in [2], df4d)
#         df4d_phase_3 = filter(row -> row.phase in [3], df4d)
#         df4d_phase_4 = filter(row -> row.phase in [4], df4d)
#         df4d_phase_5 = filter(row -> row.phase in [5], df4d)
#         df4d_phase_6 = filter(row -> row.phase in [6], df4d)
#     end
#     if fix_amp == true
#         df4d_amp_select = filter(row -> row.amp in [amplitude_select], df4d)
#     end
#     return df4d_phase_0, df4d_phase_1, df4d_phase_2, df4d_phase_3, df4d_phase_4, df4d_phase_5, df4d_phase_6, df4d_amp_select
# end

## plotting functions

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




"""
    Vars = [4, 5, 6, 9] are variables [MR,KDM5A,H4,H27,KDM6A]
"""
function Plot_C(sol)
    plt = plot(sol=sol,
        vars=[[4, 5, 6, 9]], # [MR,KDM5A,H4,H27,KDM6A]
        lw=1.5,
        xlabel="Time", ylabel="Concentration",
        # title = title,
        dpi=500)
    return plt
end





"""
    Please specifiy
    1. The Amplitude
    2. A range of frequency
    3. The db_idx of which specific dataset to calculate
    Note: the default phase ϕ = 0
"""
function anim_freq_tswitch(; range=0:0.01:0.4, amplitude=500, db_idx=301)
    t_switching_set = []
    anim = @animate for freq_i ∈ range
        _, _, _, _, ts, _, single_sol = single_solve(; db_idx=db_idx, freq=freq_i, phase=0, amplitude=amplitude, T_init=T_init, ΔT=ΔT, tspan=tspan)
        t_switching = switching_time(; sol=single_sol, pulse_period=T_init:0.1:T_init+ΔT, idx=[6, 9], return_plot=false)
        check = check_switching(single_sol, ts, tspan)
        if check == -1
            push!(t_switching_set, t_switching)
            plot(single_sol,
                # vars = [4, 5, 6, 9], # [MR,KDM5A,H4,H27]
                vars=[6, 9], # [H4,H27]
                lw=2,
                xlabel="Time", ylabel="Concentration",
                foreground_color_legend=nothing,
                # ylims = [0,500],
                title="A : $amplitude, freq = $freq_i, switch time : $t_switching",
                dpi=300)
        end
    end
    gif(anim, "freq_animation.gif", fps=1)
end




function df_freq_vs_ST_groupby_amp(df; ΔT=100, amplitude_select=false, palette=:RdYlBu_6, figure_save_path=nothing, arg...)
    @show db_idx
    if isempty(amplitude_select) == false
        df_amp_select = filter(row -> row.amp in amplitude_select, df)
        color_catg = length(amplitude_select)
    elseif isempty(amplitude_select) == true
        df_amp_select = df
    end
    
    plt = @df df_amp_select plot(
        :freq,
        :stime,
        group=:amp,
        palette=palette,
        lw=3,
        legend_title="Amplitude",
        legend_position=:topleft,
        markershape=:none,  # This removes the dots while keeping everything else the same
        dpi=500,
        arg...
    )

    xlabel!(plt, L"Driving Frequency ($\omega$)")
    ylabel!(plt, "Switching Time (ST)")
    
    if figure_save_path != nothing
        savefig(plt, figure_save_path * "ω_vs_ST" * ".png")
    end
    return plt
end










# function df4d_fix_amp_freq_vs_ST(df4d_amp_select; ΔT=100, phase_select=false, palette=:RdYlBu_6, save=false)
#     @show db_idx
#     if isempty(phase_select) == false
#         df4d_amp_select_phase_select = filter(row -> row.phase in phase_select, df4d_amp_select)
#         color_catg = length(phase_select)
#     end
#     plt = @df df4d_amp_select_phase_select plot(
#         :freq,
#         :stime,
#         group=:phase,
#         # palette = :RdYlBu_6,
#         palette=palette,
#         m=(2, 4),
#         legend_title="Phase",
#         legend_position=:outertopright,
#         # xlabel = "Switching Amplitude",
#         # ylabel = "Switching Frequency"
#         dpi=500,
#         # title = "Amplitude = 100"
#         # bg = RGB(0.2, 0.2, 0.5)
#     )
#     xlabel!(plt, "Driving Frequency")
#     ylabel!(plt, "Switching Time (ST)")
#     if save
#         class = "positive_pulse" * "_id_$db_idx" * "/ΔT = $ΔT" * "/activation_time/"
#         save_path = joinpath(@__DIR__, "figures", "switching_dynamics/$class") # generate path to save
#         isdir(save_path) || mkpath(save_path) # create directory if not exist
#         amplitude_select = unique(df4d_amp_select.amp)[1]
#         savefig(plt, save_path * "Fix_amplitude($amplitude_select)_|_ω_vs_ST=" * ".png")
#     end
#     return plt
# end



## ===== save plot and generate paths
function save_plot(plt; filename="switching_dynamics_const")
    @show db_idx
    class = "positive_pulse" * "_id_$db_idx" * "/ΔT = $ΔT" * "/activation_time/"
    # class = "oscillating_input1/" # saved two folder contains 2 ω near the transition.
    save_path = joinpath(@__DIR__, "figures", "switching_dynamics/$class") # generate path to save
    isdir(save_path) || mkpath(save_path) # create directory if not exist
    u0 = [values(u0map)[i][2] for i ∈ 1:length(u0map)]
    CSV.write(save_path * "initial_condition.csv", DataFrame(var=initi_names, values=u0))
    @show "initial_condition.csv is saved"
    CSV.write(save_path * "parameter_set.csv", DataFrame(var=vcat([p_names, "ω", "A", "ϕ"]...), values=p))
    @show "parameter_set.csv is saved"
    final_save_path = save_path * "$filename" * "_A=$amplitude" * "_w=$freq" * "_ST=$t_switching" * ".png"
    savefig(plt, final_save_path)
    @show "The saved path is here :\n $final_save_path"
end



## ====== two genes plots comparison pathgen
# function pathgen(type::String)
#     class = "2genes_example/" * type * "/"
#     save_path = joinpath(@__DIR__, "figures", "switching_dynamics/$class") # generate path to save
#     isdir(save_path) || mkpath(save_path)
#     return save_path
# end

"""
    pathgen(;num_genes::Int, type::String)
# Example 
- `type` can be `pulsatile` or `bump` or be expressive like `Dll1_vs_Dll4`
- `number_of_genes` is the number of genes
```julia
pathgen(num_genes = 2, type = "pulsatile")
pathgen(num_genes = 2, type = "Dll4_vs_Dll1")
````
"""
function pathgen(; db_idx::Union{Int,Vector{Int}}, type::String)
    # generate figure saving path
    num_genes = length(db_idx)
    identifier_1gene = "$num_genes gene/id_$db_idx/$type/"
    identifier_2genes = "$num_genes genes/id_$db_idx/$type/"
    num_genes == 1 ? identifier = identifier_1gene : identifier = identifier_2genes
    figure_save_path = joinpath(dirname(@__DIR__), "figures", identifier) # generate path to save
    isdir(figure_save_path) || mkpath(figure_save_path) #

    # generate df saving path
    df_save_path = joinpath(dirname(@__DIR__), "Data", "regular", "db_idx:$db_idx/")
    isdir(df_save_path) || mkpath(df_save_path)

    return df_save_path, figure_save_path
end

## === groupby Amplitude and frequency, and return for each (A, ω) the minimum ST 
function ST_ω(df)
    return combine(groupby(df, [:amp, :freq]), :stime => Ref => :stime, :stime => minimum)
end



# for multi gene case ---(ϕ indenpent)
function A_ω_st_multi_genes(; model=model, amplitude_range=0:50:300, freq_range=0:0.01:2, gene_set_id=[49, 592], phase_sample_size=6, prc2="NA", mute_parameter_disp=false)
    gene_set_A_ω_st_df = []
    for gene_id in gene_set_id
        gene_i_df = A_ω_st_relation(amplitude_range=amplitude_range, freq_range=freq_range, db_idx=gene_id, prc2=prc2, mute_parameter_disp=mute_parameter_disp)
        if isempty(gene_i_df)
            println("No data for gene_id = ", gene_id)
            break
        else
            @show gene_i_df
        end
        ST_ω_df = ST_ω(gene_i_df)
        ST_ω_df[!, :gene_id] .= gene_id
        @show ST_ω_df
        push!(gene_set_A_ω_st_df, ST_ω_df)
    end
    @show gene_set_A_ω_st_df
    return vcat(gene_set_A_ω_st_df...)
end


function A_ω_st_prc2_multi_genes(; gene_set_id::Vector{Int}, model=model, amplitude_range=100:20:300, freq_range=0:0.01:0.4, T_init=1e-10, ΔT=100, tspan_rt=2, prc2_range=0:0.1:1, mute_parameter_disp=true)
    all_gene_dataframes = []

    for gene_id in gene_set_id
        gene_dataframe = A_ω_st_relation_prc2_range(model=model, db_idx=gene_id, amplitude_range=amplitude_range, freq_range=freq_range, T_init=T_init, ΔT=ΔT, tspan_rt=tspan_rt, prc2_range=prc2_range, mute_parameter_disp=mute_parameter_disp)
        gene_dataframe[!, :gene_id] = fill(gene_id, nrow(gene_dataframe)) # Add a new column 'gene_id' with the current gene ID
        push!(all_gene_dataframes, gene_dataframe)
    end

    combined_gene_dataframe = vcat(all_gene_dataframes...)
    return combined_gene_dataframe
end




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


# ======== loading freq :0.0:0.02:2.0_|_amplitude :50:50:300.csv data =====
function load_csv_files(path::String, prefix::String, idx_range::AbstractRange, csv_filename::String)
    df_set = []
    for i in idx_range
        folder_path = joinpath(path, "$(prefix)$(i)")
        file_path = joinpath(folder_path, csv_filename)
        if isfile(file_path)
            df = CSV.read(file_path, DataFrame)
            push!(df_set, df)
        else
            println("File not found: $file_path")
        end
    end
    combined_df = DataFrame()
    for (i, df) in enumerate(df_set)
        df[!, :gene_id] = fill(idx_range[i], nrow(df))
        combined_df = vcat(combined_df, df)
    end
    return df_set, combined_df
end



function filter_by_gene_id(df::DataFrame, gene_id_values::Vector{Int})
    filtered_df = df[in.(df.gene_id, Ref(gene_id_values)), :]
    return filtered_df
end



function plot_prc2_vs_ST(; df_2genes::DataFrame, fixed_amp=100, freq_selection, fontsize=font(16), xlims = (0,1))
    gene_ids = unique(df_2genes.gene_id)
    fixed_amp = fixed_amp
    gene1_df = filter(row -> row.freq == freq_selection[1] && row.amp == fixed_amp && row.gene_id == gene_ids[1], df_2genes)
    gene2_df = filter(row -> row.freq == freq_selection[2] && row.amp == fixed_amp && row.gene_id == gene_ids[2], df_2genes)
    plt1 = @df gene1_df plot(:prc2, :stime,
        palette=:tab10,
        legend_title="Genes",
        legend_position=:topleft,
        linewidth=2,
        xlabel="PRC2 rate",
        ylabel="Switching Time (ST)",
        dpi=500,
        foreground_color_legend=nothing,
        title="Amplitude = $(fixed_amp)",
        label="Gene 1 (Pusatile frequency = $(freq_selection[1]))",
        tickfont=fontsize,
        guidefont=fontsize,
        legendfont=fontsize,
        xlims = xlims
    )
    plt2 = @df gene2_df plot!(plt1, :prc2, :stime,
        palette=:RdYlBu_6,
        legend_title="Genes",
        legend_position=:topleft,
        linewidth=2,
        xlabel="PRC2 rate",
        ylabel="Switching Time (ST)",
        dpi=500,
        foreground_color_legend=nothing,
        title="Amplitude = $(fixed_amp)",
        label="Gene 2 (Sustained, frequency = $(freq_selection[2]))",
        tickfont=fontsize,
        guidefont=fontsize,
        legendfont=fontsize,
        xlims = xlims
    )
end


function plot_prc2_vs_ST_amp_dependent(; df_2genes::DataFrame, fixed_amp_set=unique(df_2genes.amp), freq_selection=[0.1, 0.9], size=(1600, 1600), layout=(3, 2), fontsize=font(16), xlims = (0,1))
    plt_prc2_ST_set = []
    for fixed_amp in fixed_amp_set
        plt_prc2_ST = plot_prc2_vs_ST(df_2genes=df_2genes, fixed_amp=fixed_amp, freq_selection=freq_selection, fontsize=fontsize, xlims = xlims)
        push!(plt_prc2_ST_set, plt_prc2_ST)
    end
    return plot(plt_prc2_ST_set..., layout=layout, dpi=500, size=size)
end



#! =======================   generate Figure 9 in paper =================
function plot_Dll1vs_Dll4_prc2_ST_by_amplitude(; df_2genes=df_2genes, fixed_amp_set=[50, 300], freq_selection=[0.3, 0], layout=(1, 2), size=(800, 400), fontsize=font(8), save = false, filename = nothing)
    gr(fontfamily="Helvetica")
    Dll1vs_Dll4_prc2_ST_by_amplitude_id = plot_prc2_vs_ST_amp_dependent(; df_2genes=df_2genes, fixed_amp_set=fixed_amp_set, freq_selection=freq_selection, layout=layout, size=size, fontsize=fontsize)
    display(Dll1vs_Dll4_prc2_ST_by_amplitude_id)
    if save == true && isnothing(filename)
        savefig(Dll1vs_Dll4_prc2_ST_by_amplitude_id, "Dll1vs_Dll4_prc2_ST_by_amplitude.png")
    elseif save == true && !isnothing(filename)
        savefig(Dll1vs_Dll4_prc2_ST_by_amplitude_id, filename*".png")
    end
end