using Catalyst, DifferentialEquations, Plots, Latexify

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
ode3 = convert(ODESystem, NE_model3)
latexify(ode3)

## dynamics
#    m1      m2     k0   k1     k2     k3   d   p1  p2
p = [85.79, 1.18, 0.99, 1.71, 73.64, 19.5,  1.0, 85.8, 1.0]
u0 = rand(1:100, 5)
tspan = (0.,20.)
prob = ODEProblem(ode3,u0,tspan,p)
sol = solve(prob, Rosenbrock23())
plot(sol)

##
NE_model4 = @reaction_network begin
    (m1, m2),    M + R ↔ C
    (k0,k1),     DR0 + M + M ↔ DR1
    k2,          DR0 --> DR0 + R
    k3,          DR1 --> DR1 + R
    d,           (C, R, M) → ∅
    (p1, p2),    ∅ → (M, R)
end m1 m2 k0 k1 k2 k3 d p1 p2
@add_constraints NE_model4  begin
    DR0 + DR1  = 10
end

#ode4 = convert(ODESystem, NE_model4)
##
@show params(NE_model4)
@show species(NE_model4)
@show latexify(NE_model4)


#    m1      m2     k0   k1     k2     k3   d   p1  p2
p = [85.79, 1.18, 0.99, 1.71, 73.64, 19.5,  1.0, 85.8, 1.0]
u0 = rand(1:100, 5)
tspan = (0.,20.)
prob = ODEProblem(NE_model4,u0,tspan,p)
sol = solve(prob, Rosenbrock23())
plot(sol)








using DiffEqBiological, DifferentialEquations
##
NE_model4 = @reaction_network begin
    (m1, m2),    M + R ↔ C
    (k0,k1),     DR0 + M + M ↔ DR1
    k2,          DR0 --> DR0 + R
    k3,          DR1 --> DR1 + R
    d,           (C, R, M) → ∅
    (p1, p2),    ∅ → (M, R)
end m1 m2 k0 k1 k2 k3 d p1 p2
@add_constraints NE_model4  begin
    DR0 + DR1  = 10
end
ss = steady_states(NE_model4,p)

sort!(ss, by = x -> x[1])
sb = stability(ss,Notch_EMT,p)
sb2 = stability_tianchi(ss,Notch_EMT,p,1)






##
model = @reaction_network begin
    (n1, n2), N + R ↔ C1
    (m1, m2), M + R ↔ C2
    k0, M + M --> M + M + R
    (b0, b1), N + M ↔ C0
    d,    (N,M,R) → ∅
    p1, ∅ → N
    p2, ∅ → M
end d n1 n2 m1 m2 p1 p2 k0 b0 b1
model_p = [10, 100., 100., 50., 50., 3., 3., 20., 40., 20.,];
ss = steady_states(model, model_p)
##
using Interact
@manipulate for d = 0:0.1:100,n1= 0:0.1:100, n2= 0:0.1:100, m1= 0:0.1:100, m2= 0:0.1:100, p1= 0:0.1:100, p2= 0:0.1:100, k0= 0:0.1:100, b0= 0:0.1:100, b1= 0:0.1:100
    model_p = [d, n1, n2, m1, m2, p1, p2, k0, b0, b1]
    @show ss = steady_states(model, model_p)
end
