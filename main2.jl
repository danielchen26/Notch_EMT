## ====================== Loading packages ==========================
using DifferentialEquations
using Plots
# using ModelingToolkit
using Catalyst
using ProgressMeter
using DataFrames, CSV, Random

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
    (A * (1 + cos(w * t)), 1.0), R ↔ NR               # NICD binds RBPJ
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
end k0 k1 k2 d m p k pp kk δ α1 w A # put A at last as the control 13th variable
@show species(model)
@show parameters(model)
@show speciesmap(model)
@show paramsmap(model)





## ================================ Single Run for the model ============================================
u0 = collect(initial_condition[1, :])
p = vcat([collect(parameter_set[1, :]), 0, 0]...)
ts, cb = make_cb([50.0, 100.0], 13, 100.0)
tspan = (0.0, 150.0)
u0map = species(model) .=> u0
pmap = parameters(model) .=> p
prob = ODEProblem(model, u0map, tspan, pmap)
@time sol = solve(prob, Rosenbrock23(), callback = cb, tstops = ts)
plot(sol,
    # vars = [4, 6, 9, 5, 10], # [MR,KDM5A,H4,H27,KDM6A]
    vars = [4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
    lw = 1.5,
    xlabel = "Time", ylabel = "Concentration",
    title = "Switching Dynamics",
    dpi = 500)


##
check_switching(sol)




##  ================================ Generate A, w for swithcing events ============================================
# df = db[shuffle(1:nrow(db))[1:10], :]
parameter_set = db[:, p_names]
initial_condition = db[:, initi_names]
rand_idx_set = shuffle(1:nrow(db))
@show db_idx = rand(rand_idx_set)
@show parameter_set[db_idx, :]

##  ===========================================================================================
# Check example
# p = [5.0 ,1.0, 1.0,0.01,0.61, 6.0,15.0, 6.0,11.0, 1.0, 1.0, 2.5, 0.0]
freq = 0.0
amplitude = 300.0
p = vcat([collect(parameter_set[db_idx, :]), freq, 0.0]...)
u0 = collect(db[db_idx, initi_names])
ts, cb = make_cb([50.0, 100.0], 13, amplitude)
tspan = (0.0, 150.0)
u0map = species(model) .=> u0
pmap = parameters(model) .=> p
prob = ODEProblem(model, u0map, tspan, pmap)
@time sol = solve(prob, Rosenbrock23(), callback = cb, tstops = ts)
plot(sol,
 # vars = [4, 6, 9, 5, 10], # [MR,KDM5A,H4,H27,KDM6A]
 vars = [4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
 lw = 1.5,
 xlabel = "Time", ylabel = "Concentration",
 title = "Switching Dynamics",
 dpi = 500)
##




switch_amplitude = []
switch_frequency = []
tspan = [0.0, 150.0]

@showprogress for freq in exp10.(-4:0.1:1)
    for amplitude = 0:100
        # frequency modulation
        # p1 = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, freq, 0.0]
        # p = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, freq, 0.0]
        # u0 = [6.0, 0.0, 6.0, 40.0, 500.0, 0.6876, 0.678, 500.0, 50.6344, 1.0, 2.0]
        # u0 = [6.59, 0.0, 6.59, 43.41, 61.47, 0.02, 0.54, 48.17, 49.43, 1.09, 1.09]
        # @show parameter_set[db_idx, :]
        # @show initial_condition[db_idx, :]

        p = vcat([collect(parameter_set[db_idx, :]), freq, 0.0]...)
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
        end

        if check_switching(sol1) == -1
            break
        end
    end
end


plt = plot(switch_amplitude, switch_frequency, seriestype = :scatter, #yaxis=:log10,
    label = "switching events", title = "Frequency vs Amplitude at the switching",
    xlabel = "switch_amplitude", ylabel = "switch_frequency", dpi = 500)
add_Boundary_event(switch_amplitude, switch_frequency, plt)
# savefig(plt,"freq_vs_amplitude.png")




## ==================== Generate A, w curve for different parameter sets ============================================
rand_idx_set = shuffle(1:nrow(db))
parameter_set = db[rand_idx_set, p_names]
initial_condition = db[rand_idx_set, initi_names]


function A_w_single_run(parameter_set, initial_condition, w_range, A_range, id)
    switch_amplitude, switch_frequency = [], []
    for freq in w_range
        for amplitude in A_range
            u0 = collect(initial_condition[id, :])
            p = vcat([collect(parameter_set[id, :]), freq, 0]...)
            tspan = (0.0, 150.0)
            u0map = species(model) .=> u0
            pmap = parameters(model) .=> p
            ts, cb = make_cb([50.0, 100.0], 13, amplitude)

            prob = ODEProblem(model, u0map, tspan, pmap)
            sol = solve(prob, Rosenbrock23(), callback = cb, tstops = ts)

            check = check_switching(sol)
            if check == -1
                append!(switch_amplitude, amplitude)
                append!(switch_frequency, freq)
            end

            if check_switching(sol) == -1
                break
            end

        end
    end

    plt = plot(switch_amplitude, switch_frequency, seriestype = :scatter, #yaxis=:log10,
        label = "switching events", title = "Frequency vs Amplitude at the switching",
        xlabel = "switch_amplitude", ylabel = "switch_frequency", dpi = 500)
    return switch_amplitude, switch_frequency, plt
end
function A_w_curve(parameter_set, initial_condition, w_range, A_range, NO_sample)
    @showprogress for id = 1:NO_sample
        @show id
        switch_amplitude, switch_frequency, plt = A_w_single_run(parameter_set, initial_condition, w_range, A_range, id)
        plt_final = add_Boundary_event(switch_amplitude, switch_frequency, plt)
        save_path = joinpath(@__DIR__, "figures", "A_w_collection", "$id/") # generate path to save
        isdir(save_path) || mkpath(save_path) # create directory if not exist
        save_A_w_data(switch_amplitude, switch_frequency, save_path)
        save_init_param(parameter_set, initial_condition, id, save_path)
        # savefig(plt_final, save_path * "freq_vs_amplitude.png")
        display(plt_final)
        # display(plt)
    end
end





w_range = collect(0:0.02:3)
A_range = collect(0:1:300)
NO_sample = nrow(db)
A_w_curve(parameter_set, initial_condition, w_range, A_range, NO_sample)
