cd(@__DIR__)
using Pkg; Pkg.activate(".");
## Packages
using DiffEqFlux, DifferentialEquations, Plots; plotly()
using ModelingToolkit, DataDrivenDiffEq

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
plot(sol, vars = [H27, H0])
##

# jac = calculate_jacobian(sys)
# substitute.(jac, (u0,))
# generate_jacobian(sys, dvs = states(sys), ps = parameters(sys))
# sol[H0]



##  Model
function Notch(du,u,p,t)
  R, NR, M, MR, KDM5A, H4, H0, PRC2, H27, KDM6A, KMT = u
  k0, k1, k2, d, m, p, k, pp, kk, δ, α1, A = p
  du[1]  = dR = -A*R + NR - k1*M*R + k2*MR
  du[2]  = dNR = A*R - NR
  du[3]  = dM = -k1*M*R + k2*MR
  du[4]  = dMR = k1*M*R - k2*MR
  du[5]  = dH4 = -d*H4*KDM5A + H0*KMT
  du[6]  = dKDM5A = k0*MR + k*H27 -δ*KDM5A + α1
  du[7]  = dH0 = d*H4*KDM5A - m*H0*PRC2 + H27*KDM6A - H0*KMT
  du[8]  = dPRC2 = p*H27 - δ*PRC2 + α1
  du[9]  = dH27 = m*H0*PRC2 - H27*KDM6A
  du[10] = dKDM6A = kk*H4 - δ*KDM6A + α1
  du[11] = dKMT = pp*H4 - δ*KMT + α1
end

u0 = [6.0    ,0.0     ,6.0   ,40.0    ,500.0     ,0.6876 ,0.678  ,500.0    ,50.6344 ,1.0     ,2.0    ]
tspan = (0.0,150)
p = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, 0.0]
prob = ODEProblem(Notch, u0, tspan, p)
sol = solve(prob)
plot(sol, vars = [4,5,6,9 ], lw  = 1.5)






# ## Define DAE problem
# function Notch(out,du,u,p,t)
#   R, NR, M, MR, KDM5A, H4, H0, PRC2, H27, KDM6A, KMT = u
#   d_R, d_NR, d_M, d_MR, d_KDM5A, d_H4, d_H0, d_PRC2, d_H27, d_KDM6A, d_KMT = du
#   k0, k1, k2, d, m, p, k, pp, kk, δ, α1, A = p
#   out[1]   = -A*R + NR - k1*M*R + k2*MR - d_R
#   out[2]   = A*R - NR - d_NR
#   out[3]   = -k1*M*R + k2*MR - d_M
#   out[4]   = k1*M*R - k2*MR - d_MR
#   out[5]   = -d*H4*KDM5A + H0*KMT - d_KDM5A
#   out[6]   = k0*MR + k*H27 -δ*KDM5A + α1 - d_H4
#   out[7]   = d*H4*KDM5A - m*H0*PRC2 + H27*KDM6A - H0*KMT - d_H0
#   out[8]   = p*H27 - δ*PRC2 + α1 - d_PRC2
#   out[9]   = m*H0*PRC2 - H27*KDM6A - d_H27
#   out[10]  = kk*H4 - δ*KDM6A + α1 - d_KDM6A
#   out[11]  = pp*H4 - δ*KMT + α1  - d_KMT
#   out[12]  = M + MR - 50.0          # conservation of M
#   out[13]  = R + MR + NR - 50.0     # conservation of R
#   out[14]  = H4 + H27 + H0 - 50.0
#   nothing
# end


# u0 = [6.0, 0.0, 6.0, 40.0, 500.0, 0.6876,  0.678,  500.0,  50.6344, 1.0, 2.0, 0.0, 0.0, 0.0]
# du0 = zeros(Float64,14)
# du0[1:11] .= -0.3
# tspan = (0.0,10.0)
# differential_vars = [true,true,true,true,true,true,true,true,true,true,true,false,false,false]
# prob = DAEProblem(Notch,du0,u0,tspan,p,differential_vars=differential_vars)
# using Sundials
# sol = solve(prob)




























## SDE

# function multiplicative_noise!(du,u,p,t)
#   R, NR, M, MR, KDM5A, H4, H0, PRC2, H27, KDM6A, KMT = u
#   du[1]  = 0.0*R      #p[13]*R
#   du[2]  = 0.0*NR      #p[14]*NR
#   du[3]  = 0.0*M      #p[15]*M
#   du[4]  = 0.0*MR      #p[16]*MR
#   du[5]  = 0.0*KDM5A      #p[17]*KDM5A
#   du[6]  = 0.0*H4      #p[18]*H4
#   du[7]  = 0.0*H0      #p[19]*H0
#   du[8]  = 0.0*PRC2      #p[20]*PRC2
#   du[9]  = 0.0*H27      #p[21]*H27
#   du[10] = 0.0*KDM6A      #p[22]*KDM6A
#   du[11] = 0.0*KMT      #p[23]*KMT
# end


# p  = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, 0.0]
#
# prob = ODEProblem(Notch, u0, tspan, p)
# sol = solve(prob, Rosenbrock23())
# plot(sol)



# prob = SDEProblem(Notch,multiplicative_noise!,u0,tspan,p)
# sol = solve(prob)
# plot(sol)

# ensemble_prob = EnsembleProblem(prob)
# sol_ensemble = solve(ensemble_prob, SOSRI(), EnsembleThreads(), trajectories = 1000)
# summ  = EnsembleSummary(sol_ensemble)
# plot(summ)
