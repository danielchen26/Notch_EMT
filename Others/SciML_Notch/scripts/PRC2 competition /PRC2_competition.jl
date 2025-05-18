########################################################
using DataDrivenDiffEq
using ModelingToolkit
using DifferentialEquations
using LinearAlgebra
using Random
using Symbolics: scalarize
using Plots
using Catalyst, ProgressMeter, StatsPlots, LaTeXStrings, CSV, DataFrames


## ================= Define model  ======###########
# ===== module 1 ======== ##
function module1(; name)
    @parameters k0 k1 k2 d m p k pp kk δ α1 A t
    @variables R(t) NR(t) M(t) MR(t) KDM5A(t) H4(t) H0(t) PRC2(t) H27(t) KDM6A(t) KMT(t)
    D = Differential(t)
    sys1 = [D(R) ~ -A * R + NR - k1 * M * R + k2 * MR,
        D(NR) ~ A * R - NR,
        D(M) ~ -k1 * M * R + k2 * MR,
        D(MR) ~ k1 * M * R - k2 * MR]
    ODESystem(sys1; name)
end
@named sys_module1 = module1()

# ===== module 2 for gene 1 =========
function module2(; name)
    @parameters k0 k1 k2 d m p k pp kk δ α1 A t
    @variables R(t) NR(t) M(t) MR(t) KDM5A(t) H4(t) H0(t) PRC2(t) H27(t) KDM6A(t) KMT(t)
    D = Differential(t)
    @named common = module1()
    sys2 = [D(H4) ~ -d * H4 * KDM5A + H0 * KMT,
        D(KDM5A) ~ k0 * common.MR + k * H27 - δ * KDM5A + α1,
        D(H0) ~ d * H4 * KDM5A - m * H0 * PRC2 + H27 * KDM6A - H0 * KMT,
        D(PRC2) ~ p * H27 - δ * PRC2 + α1,
        D(H27) ~ m * H0 * PRC2 - H27 * KDM6A,
        D(KDM6A) ~ kk * H4 - δ * KDM6A + α1,
        D(KMT) ~ pp * H4 - δ * KMT + α1]
    ODESystem(sys2; name)
end
@named sys_module2_gene1 = module2()
@named sys_module2_gene2 = module2()

# ====== composing the modules ====== 
connections = [
    50 ~ sys_module2_gene1.PRC2 + sys_module2_gene2.PRC2,
    # sys_module2_gene1.KDM5A ~ -sys_module2_gene2.KDM5A + 10
]
@named connected = ODESystem(connections, t, systems=[sys_module1, sys_module2_gene1, sys_module2_gene2])
equations(connected)
u0 = [
    sys_module1.R => 0,
    sys_module1.NR => 0,
    sys_module1.M => 0,
    sys_module1.MR => 0,
    # sys_module2_gene1.common₊MR => 0,
    sys_module2_gene1.H4 => 0,
    sys_module2_gene1.KDM5A => 0,
    sys_module2_gene1.H0 => 0,
    sys_module2_gene1.PRC2 => 0,
    sys_module2_gene1.H27 => 0,
    sys_module2_gene1.KDM6A => 0,
    sys_module2_gene1.KMT => 0,
    # sys_module2_gene2.common₊MR => 0,
    sys_module2_gene2.H4 => 0,
    sys_module2_gene2.KDM5A => 0,
    sys_module2_gene2.H0 => 0,
    sys_module2_gene2.PRC2 => 0,
    sys_module2_gene2.H27 => 0,
    sys_module2_gene2.KDM6A => 0,
    sys_module2_gene2.KMT => 0
]

p = [sys_module1.A => 300,
    sys_module1.k1 => 1,
    sys_module1.k2 => 1,
    sys_module2_gene1.d => 1,
    sys_module2_gene1.k => 1,
    sys_module2_gene1.δ => 1,
    sys_module2_gene1.α1 => 1,
    sys_module2_gene1.k0 => 1,
    sys_module2_gene1.pp => 1,
    sys_module2_gene1.kk => 1,
    sys_module2_gene1.m => 1,
    sys_module2_gene1.p => 1,
    sys_module2_gene2.d => 1,
    sys_module2_gene2.k => 1,
    sys_module2_gene2.δ => 1,
    sys_module2_gene2.α1 => 1,
    sys_module2_gene2.k0 => 1,
    sys_module2_gene2.pp => 1,
    sys_module2_gene2.kk => 1,
    sys_module2_gene2.m => 1,
    sys_module2_gene2.p => 1,
]


tspan = (0.0, 100.0)
prob = ODEProblem(connected, u0, tspan, p)
sol = solve(prob, Rodas4())







connected = compose(ODESystem(connections, name=:prc2_constrained_connected), sys_module1, sys_module2_gene1, sys_module2_gene2)
connected_simp = structural_simplify(connected)
full_equations(connected_simp)
equations(connected_simp)
parameters(connected_simp)
states(connected_simp)



## import database ===========================================================
include("./Functions.jl")
db, p_names, initi_names, parameter_set, initial_condition, rand_idx_set, db_idx = loading_database(data_path="/Users/chentianchi/Desktop/M1_pro_Projects/Notch_EMT_data/Notch_params_complete.csv")

## ============================= Import old models =================================
signal_type = "bump"
model_bump = Import_model(; type=signal_type)
@show equations(model_bump)
signal_type = "pulsatile"
model_pulsatile = Import_model(; type=signal_type)
@show equations(model_pulsatile)

## =====================#### 2 genes initial consditions ###### =================================
db_idx1 = 592
freq_gene1 = 0.9
phase_gene1 = 0
p_gene1 = vcat([collect(parameter_set[db_idx1, :]), freq_gene1, 0.0, phase_gene1]...)
pmap_gene1 = parameters(model_pulsatile) .=> p_gene1
p1 = Dict(pmap_gene1)
u0_gene1 = collect(initial_condition[db_idx1, :])
u0map_gene1 = species(model_pulsatile) .=> u0_gene1
u01 = Dict(u0map_gene1)

db_idx2 = 49
freq_gene2 = 0.0
phase_gene2 = 0
p_gene2 = vcat([collect(parameter_set[db_idx2, :]), freq_gene2, 0.0, phase_gene2]...)
pmap_gene2 = parameters(model_pulsatile) .=> p_gene2
p2 = Dict(pmap_gene2)
u0_gene2 = collect(initial_condition[db_idx2, :])
u0map_gene2 = species(model_pulsatile) .=> u0_gene2
u02 = Dict(u0map_gene2)

states(connected_simp)
connected_u0map = [u01[R], u01[NR], u01[M], u01[MR], u01[H4], u01[KDM6A], u01[KMT], u02[H4], u02[KMT], u01[H0], u01[MR], u02[H0], u02[H0],
    u02[KDM5A], u02[H27], u02[MR], u02[PRC2], u02[KDM6A], u01[H0], u01[KDM6A], u02[H0], u02[KDM6A], u02[H27], u02[PRC2]]
connected_p = [p1[A], p1[k1], p1[k2], p1[d], p1[k], p1[δ], p1[α1], p1[k0], p1[m], p1[p], p1[kk], p1[pp], p2[d], p2[k], p2[δ], p2[α1], p2[k0], p2[m], p2[p], p2[kk], p2[pp]]
prob = ODAEProblem(structural_simplify(connected_simp), connected_u0map, connected_p, (0, 200.0))
sol = solve(prob, Tsit5())
plot(sol)