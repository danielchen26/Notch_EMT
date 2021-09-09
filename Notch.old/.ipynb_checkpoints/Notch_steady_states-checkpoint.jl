using DiffEqBiological, Plots, Latexify
using SymPy

# ========= Model 3 with N as the external variable
NE_model3 = @reaction_network begin
    (m1, m2),    M + R ↔ C
    (N, r),     R ↔ D
    (k0,k1),     DR0 + M + M ↔ DR1
    k2,          DR0 --> DR0 + R
    k3,          DR1 --> DR1 + R
    N,           M --> ∅
    d,           (C, D, R, M) → ∅
    (p1,p2),         ∅ → (M, R)
    end m1 m2 N r k0 k1 k2 k3 d p1 p2
# ode3 = convert(ODESystem, NE_model3) # use this in catalyst
# latexify(ode3)
# latexify(NE_model3)
@add_constraints NE_model3  begin
    DR0 + DR1  = 10
end
#   m1      m2      N      r      k0      k1      k2      k3      d      p1      p2
p = [10.39,  2.8e-2, 2.0, 1.0,  9.13e-2,  8.39,  40.47,  8.39, 1.0,  39.17,  1.0]
ss = steady_states(NE_model3,p)

@show ss




##
using Interact, DiffEqBiological
include(pwd()*"/System Biology/Julia_problems/Notch/functions.jl")

NE_model3 = @reaction_network begin
    (m1, m2),    M + R ↔ C
    (N, r),      R ↔ D
    (k0,k1),     DR0 + M + M ↔ DR1
    k2,          DR0 --> DR0 + R
    k3,          DR1 --> DR1 + R
    k4*N,        M --> ∅  # k4 = 0
    d,           (C, D, R, M) → ∅
    (p1,p2),      ∅ → (M, R)
end m1 m2 N r k0 k1 k2 k3 k4 d p1 p2

@add_constraints NE_model3  begin
    DR0 + DR1  = 1
end

# default
# model_p = [10.39, 2.8e-2, 2.0, 1.0,  9.13e-2,  8.39,  40.47,  8.39, 0.0,   1.0,  39.17,  1.0]
@manipulate for
    m1 = slider(0.01 : 0.01: 100.0, value=10.39),m2 = slider(0.01 : 0.01: 100.0, value=2.8e-2),
    N = slider(0.00 : 0.01: 100.0, value= 2.0), r = slider(0.01 : 0.01: 100.0, value=1.0),
    k0 = slider(0.01 : 0.01: 100.0, value= 0.09), k1 = slider(0.01 : 0.01: 100.0, value =8.93),
    k2 = slider(0.01 : 0.01: 100.0, value= 40.47), k3 = slider(0.01 : 0.01: 100.0, value = 8.39),
    k4 = slider(0.01 : 0.01: 100.0, value = 0.0),
    d = slider(0.01 : 0.01: 100.0, value= 1.0),
    p1 = slider(0.01 : 0.01: 100.0, value= 39.17), p2 = slider(0.01 : 0.01: 100.0, value= 1.0)

    model_p = [   m1,      m2,   N,  r,     k0,      k1,     k2,    k3,  k4,      d,    p1,   p2]

    println("params")
    df_p = DataFrame()
    df_p.variables = [:m1,:m2,:N,:r,:k0,:k1,:k2,:k3,:k4,:d,:p1,:p2];
    df_p.values = model_p
    @show df_p
    println("\n")

    ss = steady_states(NE_model3, model_p)
    df = DataFrame(ss); df.vars = NE_model3.syms
    @show df
    sort!(ss, by = x -> x[3])
    sb2 = stability_tianchi(ss,NE_model3,model_p,1)
    @show sb2
    println("\n")
end









# ## === solve steady_states
# @vars  m1 m2 N r k0 k1 k2 k3 d p1 p2 M R C D DR0 DR1
# rhs3 = [-m1*M*R + m2*C - k0*DR0*M^2 + 2*k1*DR1 - N*M - d*M + p1;
#        -m1*M*R + m2*C - N*R + r*D + k2*DR0 + k3*DR1 - d*R + p2;
#        m1*M*R - m2*C - d*C;
#        N*R - R*D - d*D
#        -k0*DR0*M^2/2 + k1*DR1;
#        k0*DR0*M^2/2 - k1*DR1];
# vars = [M, R, C, D, DR0, DR1];
#
# # # CRN rate to replace
# # p = [m1 => 10.39, m2 => 2.8e-2, N => 2.0, r => 1.0, k0 => 9.13e-2, k1 => 8.39, k2 => 40.47, k3 => 8.39, d => 1.0, p1 => 39.17, p2 => 1.0]
# # rhs3_p = [rhs3.subs(p)...]
# # S = SymPy.solve(rhs3_p, vars)


## ==========  Model 4 without N still produce bistability

NE_model4 = @reaction_network begin
    (m1, m2),    M + R ↔ C
    (k0,k1),     DR0 + M + M ↔ DR1
    k2,          DR0 --> DR0 + R
    k3,          DR1 --> DR1 + R
    d,           (C, R, M) → ∅
    (p1, p2),    ∅ → (M, R)
    end m1 m2 k0 k1 k2 k3 d p1 p2
# ode4 = convert(ODESystem, NE_model4)
# latexify(ode4)

@add_constraints NE_model4  begin
    DR0 + DR1  = 1
end
# default params
# model_p = [85.79, 1.18, 0.99, 1.71, 73.64, 19.5,  1.0, 85.8, 1.0]
@manipulate for
    m1 = slider(0.01 : 0.01: 100.0, value= 85.79), m2 = slider(0.01 : 0.01: 100.0, value = 1.18),
    k0 = slider(0.01 : 0.01: 100.0, value= 0.99), k1 = slider(0.01 : 0.01: 100.0, value= 1.71),
    k2 = slider(0.01 : 0.01: 100.0, value = 73.64), k3 = slider(0.01 : 0.01: 100.0, value =19.5),
    d = slider(0.01 : 0.01: 100.0, value = 1.0),
    p1 = slider(0.01 : 0.01: 100.0, value = 85.8), p2 = slider(0.01 : 0.01: 100.0, value = 1.0)

    model_p = [ m1, m2, k0, k1, k2, k3, d, p1, p2]

    println("params")
    df_p = DataFrame(m1  = m1, m2 = m2, k0 = k0, k1 = k1, k2 = k2, k3 = k3, d = d, p1 = p1, p2 = p2)
    @show df_p
    println("\n")

    ss = steady_states(NE_model4, model_p)
    df = DataFrame(ss); df.vars = NE_model4.syms
    @show df
    sort!(ss, by = x -> x[3])
    sb2 = stability_tianchi(ss,NE_model3,model_p,1)
    @show sb2
    println("\n")
end
