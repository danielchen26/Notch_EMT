
##
# using DifferentialEquations
using OrdinaryDiffEq
using DiffEqBiological
using DataFrames, LinearAlgebra, ProgressMeter
using Plots
using ProgressMeter
using VegaLite
using CSV

##
Notch_model_cp = @reaction_network begin
    # (N, 1.0), R ↔ NR               # NICD binds RBPJ
    # (k1, k2), M + R ↔ MR          # MITF binds RBPJ
    # k0, MR --> MR + KDM5A            # MITF-RBPJ recruit KDM5A
    d, H4 + KDM5A --> H0 + KDM5A  # Demethylation of active mark
    m, H0 + PRC2 --> H27 + PRC2   # Methylation to get repressive mark
    1.0, H27 + KDM6A --> H0 + KDM6A  # Demethylation of repressive mark
    1.0, H0 + KMT --> H4 + KMT    # Methylation to get active mark
    p, H27 --> H27 + PRC2           # PRC2 is enhenced by H27 mark 🟢
    kk, H4 --> H4 + KDM6A        # KDM6A is enhenced by H4 mark
    pp, H4 --> H4 + KMT          # KMT is enhenced by H4 mark  🟢
    k, H27 --> H27 + KDM5A                # KDM5A is enhenced by H27 mark
    δ, (PRC2, KDM5A, KDM6A, KMT) --> ∅                    # Degradation of histone reader and writers
    α1, ∅ --> (KDM6A, KMT, PRC2, KDM5A)
end d m p k pp kk δ α1 # put N at last as the control 12th variable

@add_constraints Notch_model_cp begin
    # M + MR = 50.0        # conservation of M
    # R + MR + NR = 50.0     # conservation of R
    H4 + H27 + H0 = 50.0   # conservation of H
end

##
# Define local stability function
function tianchi_stability_local(ss, model, p)

    function JE_stability(solution::Vector{Float64}, rn::DiffEqBase.AbstractReactionNetwork, p::Vector{Float64}, t = 0.0::Float64)
        jac = zeros(length(rn.syms), length(rn.syms))
        rn.jac(jac, solution, p, t)
        return (jac, eigen(jac).values)
    end
    Eigen_spectrum = [JE_stability(i, model, p)[2] for i in ss]
    tianchi_ss = [maximum(real(i)) < 1e-10 for i in Eigen_spectrum]
end

# tianchi_stability_local(ss, Notch_model_cp, model_p)

# test switch with callbacks starting from neat H27 state
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
ti = 50.0;
tf = 100.0;
t_end = 150.0;
activation = 1000.0;
ts, cb = make_cb([ti, tf], length(params(Notch_model_cp)), activation)


tt = 0
parameters_set = []
H4_H27_states = []
initial_H27_high_states = []

# # One instance of example for model_p
model_p = [0.1, 0.4, 5.0, 4.0, 4.0, 4.0, 1.0, 1.0]
steady_states(Notch_model_cp, model_p)



p_set = []
max_r = 100
@time @showprogress for d in 0.1, m in 0.4, p = 0.0:4:max_r,
    k = 0.0:4:max_r, pp = 0.0:4:max_r, kk = 0.0:4:max_r,
    δ in 1.0, α1 in 1.0

    model_p = [d, m, p, k, pp, kk, δ, α1]
    # @show model_p
    ss = steady_states(Notch_model_cp, model_p)
    sort!(ss, by = x -> x[5]) # assending H27
    # @show DataFrame(ss)
    stability = tianchi_stability_local(ss, Notch_model_cp, model_p)
    # @show stability

    # println("\n")

    if sum(stability) >= 2 && ss[end][5] > ss[end][1]# at least 2 stable fixed points
        @show DataFrame(ss)
        # ss_H27_high = ss[end]
        # @show ss_H27_high
        # start from the ss_H27_high state
        # prob = ODEProblem(Notch_model_cp, ss_H27_high, (0.0, t_end), model_p)
        # sol = solve(prob, Tsit5(), callback = cb, tstops = ts)
        # sol = solve(prob)
        # if sol(ti - 5.0)[9] > sol(ti - 5.0)[6] && sol(t_end)[6] > sol(t_end)[9]
        #     # plt = plot(sol, vars = [:MR, :KDM5A, :H4, :H27], lw = 1.5, title = "$model_p")
        #     # display(plt)
        #     # record the paramters to a DataFrame
        #     push!(parameters_set, model_p)
        #     push!(H4_H27_states, [sol(ti - 5.0)[6], sol(ti - 5.0)[9], sol(t_end)[6], sol(t_end)[9]])
        #     push!(initial_H27_high_states, ss_H27_high)
        #     push!(p_set, vcat(model_p, ["bistable", "r"]))
        # else
        push!(p_set, vcat(model_p, ["bistable", "nr"]))
        # end
    end
    push!(p_set, vcat(model_p, ["non-bistable", "NA"]))
    # if tt >= 100 # plot 10 examples
    #     break
    # end
    tt += 1
end



df = DataFrame(p_set)
df2 = DataFrame([[names(df)]; collect.(eachrow(df))], [:column; Symbol.(axes(df, 1))])[:, 2:end]
rename!(df2, [:d, :m, :p, :k, :pp, :kk, :δ, :α1, :stability, :type])
CSV.write("Notch_core_bistability_db.csv", df2)




## ==============================================================================
db = CSV.File("Notch_core_bistability_db.csv") |> DataFrame
