
## =================================================================
using ModelingToolkit
using DifferentialEquations
using Plots;
plotly();
@variables t x(t) f(t)
D = Differential(t)
@parameters A ω ϕ
f_fun(t, A, ω, ϕ) = t >= 100 && t <= 300 ? A * (1 + sign(cos(ω * t + ϕ))) : 0
@register f_fun(t, A, ω, ϕ)
# @named fol_external_f = ODESystem([f ~ A * (1 + sign(cos(ω * t))), D(D(x)) ~ -0.1*D(x) + f - 10x])
@named ldo = ODESystem([f ~ f_fun(t, A, ω, ϕ), D(D(x)) ~ -0.1D(x) + f - 10x])
sys = structural_simplify(ldo)
u0 = [D(x) => 0.0, x => 0.0]
p = [A => 3.0, ω => 0.5, ϕ => 5.5]
tspan = (0.0, 500.0)
prob = ODEProblem(sys, u0, tspan, p)
sol = solve(prob)
plot(sol, vars=[x])


pos = sol[1, :]
pos_later = pos[length(sol.t), end] # try to get the steady state solution max.
cri = maximum(abs.(pos_later)) > 1.0


