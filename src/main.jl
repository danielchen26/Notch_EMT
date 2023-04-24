
##      Environment and packages
using DifferentialEquations
using Plots
using ModelingToolkit
using Catalyst
using ProgressMeter
##

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


## =================== Model =================================
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

# initial condition [R(t), NR(t), M(t), MR(t), KDM5A(t), H4(t), H0(t), PRC2(t), H27(t), KDM6A(t), KMT(t)]
u0 = [6.59, 0.0, 6.59, 43.41, 61.47, 0.02, 0.54, 48.17, 49.43, 1.09, 1.09]


u0 = [6.0, 0.0, 6.0, 40.0, 500.0, 0.6876, 0.678, 500.0, 50.6344, 1.0, 2.0]


# # with w =0, when N is turned on, it is a constant rate
p = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, 0.0, 0.0]

# with w = 1,when N is turned on, it is a oscillating rate
p = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, 0.05, 0.0]

ts, cb = make_cb([50.0, 100.0], 13, 180.0)
# prob = ODEProblem(model, u0, tspan, p)


u0 = [6.6, 0.0, 6.6, 43.4, 0.2, 948.3, 0.7, 1997.3, 399.4, 1.7, 1.2]
p = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, 1.0, 0.0]


# mappings from symbolic variables to their values
u0map = species(model) .=> u0
pmap = parameters(model) .=> p
tspan = [0.0, 150.0]
prob = ODEProblem(model, u0map, tspan, pmap)
@time sol = solve(prob, Rosenbrock23(), callback = cb, tstops = ts)
plot(sol,
    # vars = [4, 6, 9, 5, 10], # [MR,KDM5A,H4,H27,KDM6A]
    vars = [6, 9], # [MR,KDM5A,H4,H27,KDM6A]
    lw = 1.5,
    xlabel = "Time", ylabel = "Concentration",
    title = "Switching Dynamics")

## too slow

# u0 = [6.0    ,0.0     ,6.0   ,40.0    ,500.0     ,0.6876 ,0.678  ,500.0    ,50.6344 ,1.0     ,2.0    ]
# tspan = [0.0, 150.0]
# @manipulate for freq =  slider(0.00 : 0.01: 10.0, value= 1.0,   label ="Notch frequency"),
# 	Amplitude = slider(0.00 : 10: 1000.0, value= 1000.0,   label ="Notch Amplitude")
#
# 	p = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, freq, 0.0]
# 	ts, cb  = make_cb([50.0,100.0], 13, Amplitude)
# 	prob = ODEProblem(model, u0, tspan, p)
# 	sol = solve(prob, Rosenbrock23(), callback=cb, tstops=ts)
# 	plot(sol, vars = [5,6,9,10], lw  = 1.5)
#
# end

##

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
check_switching(sol)


switch_amplitude = []
switch_frequency = []
@showprogress for freq = 0:0.1:2, amplitude = 0:50
    # frequency modulation
    p = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, freq, 0.0]
    pmap = parameters(model) .=> p
    u0 = [6.0, 0.0, 6.0, 40.0, 500.0, 0.6876, 0.678, 500.0, 50.6344, 1.0, 2.0]
    u0map = species(model) .=> u0
    ts, cb = make_cb([50.0, 100.0], 13, amplitude)
    prob1 = ODEProblem(model, u0map, tspan, pmap)
    sol1 = solve(prob1, Rosenbrock23(), callback = cb, tstops = ts)
    plt1 = plot(sol1, vars = [5, 6, 9, 10], lw = 1.5, title = "Amp: $amplitude, Frequency: $freq")
    display(plt1)
    check = check_switching(sol1)
    if check == -1
        append!(switch_amplitude, amplitude)
        append!(switch_frequency, freq)
    end
end


plt = plot(switch_amplitude, switch_frequency, seriestype = :scatter, #yaxis=:log10,
    label = "switching events", title = "Frequency vs Amplitude at the switching",
    xlabel = "switch_amplitude", ylabel = "switch_frequency", dpi = 500)
# savefig(plt,"freq_vs_amplitude.png")



##
using StatsPlots
@df df marginalscatter(:switch_amplitude, :switch_frequency)



using DataFrames, CSV
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

add_Boundary_event(switch_amplitude, switch_frequency, plt)



##
for freq = 0.05:0.1:1.05, amplitude = 30:20:190
    # frequency modulation
    p1 = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, freq, 0.0]
    ts, cb = make_cb([50.0, 100.0], 13, amplitude)
    prob1 = ODEProblem(model, u0, tspan, p1)
    sol1 = solve(prob1, Rosenbrock23(), callback = cb, tstops = ts)
    plt1 = plot(sol1, vars = [5, 6, 9, 10], lw = 1.5, title = "Amp: $amplitude, Frequency: $freq")

    # amplitude modulation
    p2 = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, 0.0, 0.0]
    ts, cb = make_cb([50.0, 100.0], 13, amplitude)
    prob2 = ODEProblem(model, u0, tspan, p2)
    sol2 = solve(prob2, Rosenbrock23(), callback = cb, tstops = ts)
    plt2 = plot(sol2, vars = [5, 6, 9, 10], lw = 1.5, title = "Amp: $amplitude")


    plt = plot(plt1, plt2)
    display(plt)
    # savefig(plt, "./figures/freq = 0.05:0.1:2.05, amplitude = 30:20:190/frequency_$freq" * " Amplitude_$amplitude.png")
end

##
# freq = 3;
# amplitude = 50;
# p1 = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, freq, 0.0]
# ts, cb = make_cb([50.0, 100.0], 13, amplitude)
# prob1 = ODEProblem(model, u0, tspan, p1)
# sol1 = solve(prob1, Rosenbrock23(), callback = cb, tstops = ts)
# plt1 = plot(sol1, vars = [2, 5, 6, 9, 10], lw = 1.5, title = "Amp: $amplitude, Frequency: $freq")


# p2 = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, 0.0, 0.0]

# ts, cb = make_cb([50.0, 100.0], 13, amplitude)
# prob2 = ODEProblem(model, u0, tspan, p2)
# sol2 = solve(prob2, Rosenbrock23(), callback = cb, tstops = ts)
# plt2 = plot(sol2, vars = [5, 6, 9, 10], lw = 1.5, title = "Amp: $amplitude")

# plt = plot(plt1, plt2, size = (800, 400), dpi = 500)
# display(plt)
# savefig(plt, "./figures/frequency_$freq" * " Amplitude_$amplitude.png")

##
# dprob = DiscreteProblem(model, u0, tspan, p1)
# jprob = JumpProblem(model, dprob, Direct())
# jsol = solve(jprob, callback=cb, tstops=ts,SSAStepper())
# jplt = plot(jsol, vars = [ :KDM5A,:MR, :H4, :H27], lw  = 1.5)
#
# plot(dplt,jplt, layout = (2,1), figsize = (800,600))

## animation
for freq ∈ 0.05:0.1:1.05
    anim = @animate for amplitude ∈ 0:10:200

        p1 = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, freq, 0.0]
        ts, cb = make_cb([50.0, 100.0], 13, amplitude)
        prob1 = ODEProblem(model, u0, tspan, p1)
        sol1 = solve(prob1, Rosenbrock23(), callback = cb, tstops = ts)
        plt1 = plot(sol1, vars = [5, 6, 9, 10], lw = 1.5, title = "Amp: $amplitude, Frequency: $freq")

        p2 = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, 0.0, 0.0]

        ts, cb = make_cb([50.0, 100.0], 13, amplitude)
        prob2 = ODEProblem(model, u0, tspan, p2)
        sol2 = solve(prob2, Rosenbrock23(), callback = cb, tstops = ts)
        plt2 = plot(sol2, vars = [5, 6, 9, 10], lw = 1.5, title = "Amp: $amplitude")

        plt = plot(plt1, plt2, size = (800, 400), dpi = 500)
        display(plt)
    end
    gif(anim, "figures/gif/freq=$freq.gif", fps = 1)
end






## exploring phase dependence of epigenetic switching at the critical_freq
