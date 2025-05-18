# cd(@__DIR__)
# using Pkg; Pkg.activate(".")
## Generate some data by solving a differential equation
########################################################
using DataDrivenDiffEq
using ModelingToolkit
using DifferentialEquations
using LinearAlgebra
using Random
using Symbolics: scalarize
using Plots


## ================= Define model  ======###########
@parameters k0 k1 k2 d m p k pp kk δ α1 A t
@variables R(t) NR(t) M(t) MR(t) KDM5A(t) H4(t) H0(t) PRC2(t) H27(t) KDM6A(t) KMT(t)
D = Differential(t)


##
eqs = [ D(R)     ~ -A*R + NR - k1*M*R + k2*MR,
        D(NR)    ~ A*R - NR,
        D(M)     ~ -k1*M*R + k2*MR,
        D(MR)    ~ k1*M*R - k2*MR,
        D(H4)    ~ -d*H4*KDM5A + H0*KMT,
        D(KDM5A) ~ k0*MR + k*H27 - δ*KDM5A + α1,
        D(H0)    ~ d*H4*KDM5A - m*H0*PRC2 + H27*KDM6A - H0*KMT,
        D(PRC2)  ~ p*H27 - δ*PRC2 + α1,
        D(H27)   ~ m*H0*PRC2 - H27*KDM6A,
        D(KDM6A) ~ kk*H4 - δ*KDM6A + α1,
        D(KMT)   ~ pp*H4 - δ*KMT + α1  ]
@named sys = ODESystem(eqs)
##
u0 = [R =>6.0, NR => 0.0, M => 6.0, MR => 40.0, KDM5A => 500.0,  H4 =>0.6876,
      H0 => 0.678, PRC2 => 500.0,  H27 => 50.6344 , KDM6A =>1.0, KMT =>2.0 ]
model_p = [k0 => 10.0, k1 =>1.0, k2 =>1.0, d =>0.2, m => 0.53,
           p => 1.8, k => 3.77, pp => 19.08, kk => 19.08, δ => 1.0, α1 => 1.0, A =>0.0]
tspan = (0.0,100.0)

prob = ODEProblem(sys,u0,tspan,model_p,jac=true)
sol = solve(prob,Tsit5())



## Start the automatic discovery
ddprob = ContinuousDataDrivenProblem(sol)

# @variables t R(t) NR(t) M(t) MR(t) KDM5A(t) H4(t) H0(t) PRC2(t) H27(t) KDM6A(t) KMT(t)
u = [R;NR;M;MR;KDM5A;H4;H0;PRC2;H27;KDM6A;KMT]
basis = Basis(polynomial_basis(u, 2), u, iv = t)
opt = STLSQ(exp10.(-5:0.1:-1))
ddsol = solve(ddprob, basis, opt, normalize = true, progress = true) # did not finish
print(ddsol, Val{true})




## === module 1 ====
module1 = [ D(R)     ~ -A*R + NR - k1*M*R + k2*MR,
            D(NR)    ~ A*R - NR,
            D(M)     ~ -k1*M*R + k2*MR,
            D(MR)    ~ k1*M*R - k2*MR]
@named sys_m1 = ODESystem(module1)
##
u0_m1 = [R =>6.0, NR => 0.0, M => 6.0, MR => 40.0]
m1_p = [k1 =>1.0, k2 =>1.0, A =>0.0]
tspan = (0.0,100.0)
prob = ODEProblem(sys_m1,u0_m1,tspan,m1_p,jac=true)
sol = solve(prob,Tsit5())

## Start the automatic discovery
ddprob = ContinuousDataDrivenProblem(sol)
u = [R;NR;M;MR]
w = [k1;k2;A]
basis = Basis(polynomial_basis(u, 5), u, parameters = w, iv = t)
opt = STLSQ(exp10.(-5:0.1:-1))
@time ddsol = solve(ddprob, basis, opt, normalize = true, progress =  true)
print(ddsol, Val{true})




## ====== example 1 ===== #
Random.seed!(1111) # For the noise

# Create a nonlinear pendulum
function pendulum(u, p, t)
    x = u[2]
    y = -9.81sin(u[1]) - 0.3u[2]^3 -3.0*cos(u[1]) - 10.0*exp(-((t-5.0)/5.0)^2)
    return [x;y]
end

u0 = [0.99π; -1.0]
tspan = (0.0, 15.0)
prob = ODEProblem(pendulum, u0, tspan)
sol = solve(prob, Tsit5(), saveat = 0.01)

# Create the data with additional noise
X = sol[:,:] + 0.1 .* randn(size(sol))
DX = similar(sol[:,:])

for (i, xi) in enumerate(eachcol(sol[:,:]))
    DX[:,i] = pendulum(xi, [], sol.t[i])
end

ts = sol.t
##
prob = ContinuousDataDrivenProblem(X, ts, GaussianKernel() ,
    U = (u,p,t)->[exp(-((t-5.0)/5.0)^2)], p = ones(2))

@variables u[1:2] c[1:1]
@parameters w[1:2]
u = scalarize(u)
c = scalarize(c)
w = scalarize(w)

h = Num[sin.(w[1].*u[1]);cos.(w[2].*u[1]); polynomial_basis(u, 5); c]

basis = Basis(h, u, parameters = w, controls = c)
##
λs = exp10.(-10:0.1:-1)
opt = STLSQ(λs)
res = solve(prob, basis, opt, progress = true, denoise = false, normalize = false, maxiter = 500)
system = result(res)
println(system)
params = parameters(res)
ps = parameter_map(res)
println(ps)










## === example 2 =================================
Random.seed!(1112) # For the noise

# Create a nonlinear pendulum
# function notch!(du, u, p, t)
#         R, NR, M, MR, KDM5A, H4, H0, PRC2, H27, KDM6A, KMT = u
#         k0, k1, k2, d, m, p, k, pp, kk, δ, α1, A  = p
#         du[1]  = d_R     =   -A*R + NR - k1*M*R + k2*MR
#         du[2]  = d_NR    =   A*R - NR
#         du[3]  = d_M     =   -k1*M*R + k2*MR
#         du[4]  = d_MR    =   k1*M*R - k2*MR
#         du[5]  = d_H4    =   -d*H4*KDM5A + H0*KMT
#         du[6]  = d_KDM5A =   k0*MR + k*H27 - δ*KDM5A + α1
#         du[7]  = d_H0    =   d*H4*KDM5A - m*H0*PRC2 + H27*KDM6A - H0*KMT
#         du[8]  = d_PRC2  =   p*H27 - δ*PRC2 + α1
#         du[9]  = d_H27   =   m*H0*PRC2 - H27*KDM6A
#         du[10] = d_KDM6A =   kk*H4 - δ*KDM6A + α1
#         du[11] = d_KMT   =   pp*H4 - δ*KMT + α1
# end

# u0 = [6.0,  0.0,  6.0,  40.0,  500.0,  0.6876,  0.678,  500.0,   50.6344 , 1.0, 2.0 ]
# p_ = [ 10.0, 1.0, 1.0, 0.2,  0.53,  1.8,  3.77,  19.08,  19.08,  1.0,  1.0, 0.0]
# tspan = (0.0, 15.0)
# prob = ODEProblem(notch!, u0, tspan, p_)
# sol = solve(prob)
# plot(sol)

@parameters k0 k1 k2 d m p k pp kk δ α1 A t
@variables R(t) NR(t) M(t) MR(t) KDM5A(t) H4(t) H0(t) PRC2(t) H27(t) KDM6A(t) KMT(t)
D = Differential(t)


##
eqs = [ D(R)     ~ -A*R + NR - k1*M*R + k2*MR,
        D(NR)    ~ A*R - NR,
        D(M)     ~ -k1*M*R + k2*MR,
        D(MR)    ~ k1*M*R - k2*MR,
        D(H4)    ~ -d*H4*KDM5A + H0*KMT,
        D(KDM5A) ~ k0*MR + k*H27 - δ*KDM5A + α1,
        D(H0)    ~ d*H4*KDM5A - m*H0*PRC2 + H27*KDM6A - H0*KMT,
        D(PRC2)  ~ p*H27 - δ*PRC2 + α1,
        D(H27)   ~ m*H0*PRC2 - H27*KDM6A,
        D(KDM6A) ~ kk*H4 - δ*KDM6A + α1,
        D(KMT)   ~ pp*H4 - δ*KMT + α1  ]
@named sys = ODESystem(eqs)
##
u0 = [R =>6.0, NR => 0.0, M => 6.0, MR => 40.0, KDM5A => 500.0,  H4 =>0.6876,
      H0 => 0.678, PRC2 => 500.0,  H27 => 50.6344 , KDM6A =>1.0, KMT =>2.0 ]
model_p = [k0 => 10.0, k1 =>1.0, k2 =>1.0, d =>0.2, m => 0.53,
           p => 1.8, k => 3.77, pp => 19.08, kk => 19.08, δ => 1.0, α1 => 1.0, A =>0.0]
tspan = (0.0,100.0)

prob = ODEProblem(sys,u0,tspan,model_p,jac=true)
sol = solve(prob,Tsit5())

plot(sol)

##
# Create the data with additional noise
X = sol[:,:] + 0.1 .* randn(size(sol))

ts = sol.t

##
ddprob = ContinuousDataDrivenProblem(X, ts, GaussianKernel())

u = [R;NR;M;MR;KDM5A;H4;H0;PRC2;H27;KDM6A;KMT]
Ψ = Basis(polynomial_basis(u, 2), u, iv = t)
opt = STLSQ(exp10.(-10:0.1:1))
res = solve(ddprob, Ψ, opt, normalize = true, progress =  true)
print(res, Val{true})

system = result(res)
params = parameters(res)

infered_prob = ODEProblem(system, u0, tspan, parameters(res))
infered_solution = solve(infered_prob, Tsit5(), saveat = ts)





🍏 ## working on this right now ============== subsystem =============================================================================
sub2 = [D(H4)    ~ -d*H4*KDM5A + H0*KMT,
        D(KDM5A) ~ k0 + k*H27 - δ*KDM5A + α1,
        D(H0)    ~ d*H4*KDM5A - m*H0*PRC2 + H27*KDM6A - H0*KMT,
        D(PRC2)  ~ p*H27 - δ*PRC2 + α1,
        D(H27)   ~ m*H0*PRC2 - H27*KDM6A,
        D(KDM6A) ~ kk*H4 - δ*KDM6A + α1,
        D(KMT)   ~ pp*H4 - δ*KMT + α1 ]
@named sys_sub2 = ODESystem(sub2)

u0 = [KDM5A => 500.0,  H4 =>0.6876,
      H0 => 0.678, PRC2 => 500.0,  H27 => 50.6344 , KDM6A =>1.0, KMT =>2.0 ]
model_p = [k0 => 10.0, k1 =>1.0, k2 =>1.0, d =>0.2, m => 0.53,
           p => 1.8, k => 3.77, pp => 19.08, kk => 19.08, δ => 1.0, α1 => 1.0, A =>0.0]
tspan = (0.0,30.0)

prob = ODEProblem(sys_sub2,u0,tspan,model_p,jac=true)
sol = solve(prob, Tsit5(), saveat=0.005)

plot(sol)
##
# Create the data with additional noise
X = sol[:,:] + 0.1 .* randn(size(sol))

ts = sol.t

##
ddprob = ContinuousDataDrivenProblem(X, ts, GaussianKernel())
ddprob = ContinuousDataDrivenProblem(sol)



u = [KDM5A;H4;H0;PRC2;H27;KDM6A;KMT]
w = [k0 k1 k2 d m p k pp kk δ α1 A t]
Ψ = Basis(polynomial_basis(u, 2), u, iv = t)
Ψ = Basis([H4*KDM5A,H0*KMT,H0*PRC2,H27*KDM6A,H0*KMT,H0*PRC2,H27*KDM6A], [u], iv = t)
Ψ = Basis([u;H4*KDM5A;H0*KMT;H0*PRC2;H27*KDM6A;H0*KMT;H0*PRC2;H27*KDM6A], u, iv = t)
opt = STLSQ(exp10.(-7:0.1:-1))
res = solve(ddprob, Ψ, opt, normalize = true, progress =  true)
print(res, Val{true})

system = result(res)
params = parameters(res)


u00 = [ 500.0,  0.6876, 0.678,  500.0,   50.6344 , 1.0, 2.0 ]
infered_prob = ODEProblem(system, u00, tspan, parameters(res))
infered_solution = solve(infered_prob, Tsit5())

plot(infered_solution)
