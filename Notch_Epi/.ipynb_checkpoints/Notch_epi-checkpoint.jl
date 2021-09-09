# -*- coding: utf-8 -*-
# ##
# using ModelingToolkit, OrdinaryDiffEq
# using DiffEqBiological, LinearAlgebra
# using Plots, Latexify
#
#
# ##
# notch = @reaction_network begin
#     (δ, k0),           N ↔ ∅
#     (k1, k2),          N + R ↔ NR
#     (kf, kb),          R + D ↔ RD
#      γ,                HW + D → HW + DA
#  end δ k0 k1 k2 kf kb γ
# latexify(notch)
#
# p = [1, 0.2, 1., 0.25, 0.25, 0.22, 2]
# u0=[]
#
# @add_constraints notch begin
#     D + DA + RD= 10
# end
# ss = steady_states(notch,p)
# stability(ss,notch_emt,p)


##
using Interact, DiffEqBiological, DataFrames
include(pwd()*"/Notch_Epi/functions.jl")

Notch_model = @reaction_network begin

# 	The competitve binding between NICD and MITF
# 	model NICD as an external input

  (N,k0),              R ↔ NR               # NICD binds RBPJ
  (k1, k2),            M + R ↔ MR 		    # MITF binds RBPJ
   # hill(M,v1,Kd1,2),   M --> M + R          # MITF form homodimer to activate RBPJ
   # hillr(N,v2,Kd2,1),   N + M --> N         # NICD transcriptionally inhibit MITF

# 	The epigenetic core regulation for the histone state of mir-222

	d1, 		H4  + KDM5A  --> H0  + KDM5A  # Demethylation of active mark
	m0,			H0  + PRC2   --> H27 + PRC2   # Methylation to get repressive mark
	d0, 		H27 + KDM6A  --> H0  + KDM6A  # Demethylation of repressive mark
	m1, 		H0  + KMT    --> H4  + KMT    # Methylation to get active mark

# 	Epigenetic feeback
	k0,			MR --> MR + KDM5A    				# MITF-RBPJ recruit KDM5A
	H27, 		∅ --> PRC2           				# PRC2 is enhenced by H27 mark
	H4,	 		∅ --> KDM6A        					# KDM6A is enhenced by H4 mark
	H4, 		∅ --> KMT          					# KMT is enhenced by H4 mark
	H27, 		∅ --> KDM5A        				    # KDM5A is enhenced by H27 mark
	δ,			(PRC2, KDM5A, KDM6A, KMT) --> ∅     # Degradation of histone reader and writers
# basal creation
	α1,          ∅ --> (KDM5A, KDM6A, KMT)
	α2,          ∅ --> PRC2
end N k0 k1 k2 d1 d0 m1 m0 δ α1 α2 #v1 Kd1 v2 Kd2
@add_constraints Notch_model  begin
    R + NR + MR  = 10.0
	# H4 + H0 + H27 = 100.0
end

# default
# model_p = [10.39, 2.8e-2, 2.0, 1.0,  9.13e-2,  8.39,  40.47,  8.39, 0.0,   1.0,  39.17,  1.0]
# model_p = [N, k0, k1, k2, d1, d0, m1, m0, δ, v1, Kd1, v2, Kd2]
model_p = [0.0, 10.0, 10.0, 10.0, 100.0, 100.0, 100.0, 100.0, 1.0, 1.0, 50.0, 1.0, 50.0, 1.0, 1.9]
steady_states(Notch_model, model_p)


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








# #################################################################
# Use the following
#  #################################################################

## =================================================================
using Interact, DiffEqBiological, DataFrames
include(pwd()*"/Notch_Epi/functions.jl")
Notch_model2 = @reaction_network begin
  (k1, k2),            M + R ↔ MR 		    # MITF binds RBPJ
	d1, 		H4  + KDM5A  --> H0  + KDM5A  # Demethylation of active mark
	m0,			H0  + PRC2   --> H27 + PRC2   # Methylation to get repressive mark
	d0, 		H27 + KDM6A  --> H0  + KDM6A  # Demethylation of repressive mark
	m1, 		H0  + KMT    --> H4  + KMT    # Methylation to get active mark

	k0,			MR --> MR + KDM5A    				        # MITF-RBPJ recruit KDM5A
	p, 		    H27 --> H27 + PRC2           				# PRC2 is enhenced by H27 mark
	kk,	 		H4 --> H4 + KDM6A        					# KDM6A is enhenced by H4 mark
	pp, 		H4 --> H4 + KMT          					# KMT is enhenced by H4 mark
	k, 		    H27 --> H27 + KDM5A        				    # KDM5A is enhenced by H27 mark
	δ,			(PRC2, KDM5A, KDM6A, KMT) --> ∅             # Degradation of histone reader and writers

	α1,          ∅ --> (KDM5A, KDM6A, KMT)
	α2,          ∅ --> PRC2
end k0 k1 k2 d1 d0 m1 m0 p k pp kk δ α1 α2

@add_constraints Notch_model2  begin
    M + MR  = 10.0
	R + MR  = 10.0
	H4 + H27 + H0 = 10.0
end

##
@show params(Notch_model2)
model_p = [106.8, 137.62, 11.73, 0.913, 60.215, 27.57, 3.927, 223.89, 223.89, 27.57, 27.57, 1.0, 1.0, 1.98]
ss = steady_states(Notch_model2, model_p)
##
@manipulate for
    k0 = slider(0.00 : 0.01: 100.0, value=106.8),
    k1 = slider(0.00 : 0.01: 100.0, value= 137.62),
    k2 = slider(0.00 : 0.01: 100.0, value= 11.73),
    d1 = slider(0.00 : 0.01: 100.0, value= 0.913),
    d0 = slider(0.00 : 0.01: 100.0, value = 60.215),
    m1 = slider(0.00 : 0.01: 100.0, value= 27.57),
    m0 = slider(0.00 : 0.01: 100.0, value= 3.927),
	p = slider(0.00 : 0.01: 100.0, value= 223.89),
	k = slider(0.00 : 0.01: 100.0, value=223.89),
	pp = slider(0.0 : 0.01: 100.0, value=27.57),
	kk = slider(0.0 : 0.01: 100.0, value = 27.57),
 	δ = slider(0.0 : 0.01: 100.0, value = 1.00),
 	α1 =  slider(0.0 : 0.01: 100.0, value= 1.00),
	α2 = slider(0.0 : 0.01: 100.0, value= 1.98)

    model_p = [ k0,  k1,  k2,  d1,  d0,  m1,  m0,  p,  k,  pp,  kk,  δ,  α1,  α2];

    # println("params")
    # df_p = DataFrame();
    # df_p.variables = [:M,:R,:MR,:H4,:KDM5A,:H0,:PRC2,:H27,:KDM6A,:KMT];
    # df_p.values = model_p;
    # @show df_p
    # println("\n")

    ss = steady_states(Notch_model2, model_p)
    df = DataFrame(ss); df.vars = Notch_model2.syms
    @show df
    sort!(ss, by = x -> x[4])
    sb2 = stability_tianchi(ss,Notch_model2,model_p,1)
    # sb2 = stability(ss,Notch_model2,model_p)
    @show sb2
    println("\n")
end

## ODE simulation
using DifferentialEquations,Plots;plotly()
Notch_rd_ODE = @ode_def begin
	dM =  k2*MR - M*R*k1
	dR =   k2*MR - M*R*k1
	dMR =  -k2*MR + M*R*k1
	dH4 =  -H4*d1*KDM5A + m1*H0*KMT
	dKDM5A =  α1 + k*H27 + k0*MR - δ*KDM5A
	dH0 = H4*d1*KDM5A + d0*H27*KDM6A - m0*H0*PRC2 - m1*H0*KMT
	dPRC2 =   α2 + p*H27 - δ*PRC2
	dH27 =  -d0*H27*KDM6A + m0*H0*PRC2
	dKDM6A =  α1 + kk*H4 - δ*KDM6A
	dKMT =   α1 + pp*H4 - δ*KMT
end k0 k1 k2 d1 d0 m1 m0 p k pp kk δ α1 α2÷
model_p = [106.8, 137.62, 11.73, 0.913, 60.215, 27.57, 3.927, 223.89, 223.89, 27.57, 27.57, 1.0, 1.0, 1.98]
u0 = rand(0.0:100.0, size(Notch_rd_ODE.syms))
tspan = (0.0,30)
prob = ODEProblem(Notch_rd_ODE,u0,tspan,model_p)
sol = solve(prob)
plot(sol, vars = [:MR, :H4, :H27])

##
function initialization_con(H_c, M_c, R_c)
	# inital conditions
	H_set = randfixsum(1, 3, H_c)
	M_con = M_c; R_con = R_c;
	MR = rand(0:minimum([M_c, R_c]),1)[1]
	M = M_con - MR
	R = R_con -MR
	H4, H0, H27 = H_set[1][1], H_set[1][2], H_set[1][3]
	KDM5A, PRC2, KDM6A, KMT = rand(0:100,4)
	u0 = [M,R,MR,H4,KDM5A,H0,PRC2,H27,KDM6A,KMT]
end

plt =plot()
for inital = 1:100
	# u0 = rand(0.0:100.0, size(Notch_rd_ODE.syms))
	u0 = initialization_con(10.0, 10.0, 10.0)
	tspan = (0.0,1e2)
	prob = ODEProblem(Notch_rd_ODE,u0,tspan,model_p)
	sol = solve(prob)
	plot!(plt, sol, vars = [:MR, :H4, :H27], linecolor = [:orange :green :red], legend = false)
	# plot!(plt, sol, vars = [:H4],  legend = false, lw =1 )
end

## calculate the MR, H4, H27 equilibrium points
result =[]
u0_set = []
for repeat = 1:10000
	# inital conditions
	u0 = initialization_con(50.0, 50.0, 50.0)
	push!(u0_set,u0)
	tspan = [0,1e2]
	prob = ODEProblem(Notch_rd_ODE,u0,tspan,model_p)
	sol = solve(prob)
	push!(result, sol[end][[3,4,8]])
	# push!(result, sol[end])
end
rr = [round.(i,digits = 4) for i in result]
state_var = DataFrame(hcat(rr...)'); rename!(state_var,[:MR, :H4, :H27])
u0_db = DataFrame(hcat(u0_set...)'); rename!(u0_db,[:u0_M, :u0_R, :u0_MR, :u0_H4, :u0_KDM5A, :u0_H0, :u0_PRC2, :u0_H27, :u0_KDM6A, :u0_KMT])
db = hcat([state_var,u0_db]...)
unique(rr)



## #################################
# The complete model
# I only want to  model the parameters of the relative methylation
####################################
using Interact, DiffEqBiological, DataFrames
# include(pwd()*"/Notch_Epi/functions.jl")
include(pwd()*"/functions.jl")
Notch_model_cp = @reaction_network begin
    (N,1.0),              R ↔ NR               # NICD binds RBPJ
    (k1, k2),             M + R ↔ MR 		    # MITF binds RBPJ
     k0,		          MR --> MR + KDM5A    				        # MITF-RBPJ recruit KDM5A

	d, 	      	H4  + KDM5A  --> H0  + KDM5A  # Demethylation of active mark
	m,			H0  + PRC2   --> H27 + PRC2   # Methylation to get repressive mark
	1.0, 		H27 + KDM6A  --> H0  + KDM6A  # Demethylation of repressive mark
	1.0, 		H0  + KMT    --> H4  + KMT    # Methylation to get active mark

	p, 		    H27 --> H27 + PRC2           				# PRC2 is enhenced by H27 mark
	kk,	 		H4 --> H4 + KDM6A        					# KDM6A is enhenced by H4 mark
	pp, 		H4 --> H4 + KMT          					# KMT is enhenced by H4 mark
	k, 		    H27 --> H27 + KDM5A        				    # KDM5A is enhenced by H27 mark
	δ,			(PRC2, KDM5A, KDM6A, KMT) --> ∅             # Degradation of histone reader and writers

	α,          ∅ --> (KDM5A, KDM6A, KMT, PRC2)
end k0 k1 k2 d m p k pp kk δ α N

@add_constraints Notch_model_cp  begin
    M + MR  = 500.0        # conservation of M
	R + MR + NR = 50.0     # conservation of R
	H4 + H27 + H0 = 50.0   # conservation of H
end

@show params(Notch_model_cp)
model_p = [106.8, 137.62, 11.73, 0.913, 3.927, 223.89, 223.89, 27.57, 27.57, 1.0, 1.0, 0.0]
ss = steady_states(Notch_model_cp, model_p)


##
@manipulate for
    k0 = slider(0.00 : 0.01: 100.0, value=106.8),     # MR --> KDM5A
    k1 = slider(0.00 : 0.01: 100.0, value= 137.62),   # M + R ↔ MR forward
    k2 = slider(0.00 : 0.01: 100.0, value= 11.73),    # M + R ↔ MR backward
    d = slider(0.00 : 0.01: 100.0, value= 0.913),     # Demethylation of active mark
    m = slider(0.00 : 0.01: 100.0, value= 3.927),     # Methylation to get repressive mark
	p = slider(0.00 : 0.01: 100.0, value= 223.89),    # PRC2 is enhenced by H27 mark
	k = slider(0.00 : 0.01: 100.0, value=223.89),     # KDM5A is enhenced by H27 mark
	pp = slider(0.0 : 0.01: 100.0, value=27.57),      # KMT is enhenced by H4 mark
	kk = slider(0.0 : 0.01: 100.0, value = 27.57),    # KDM6A is enhenced by H4 mark
 	δ = slider(0.0 : 0.01: 100.0, value = 1.00),      # Degradation of histone reader and writers
 	α =  slider(0.0 : 0.01: 100.0, value= 1.00),      # Production rate
    N =  slider(0.0 : 0.01: 100.0, value= 0.00)       # Production rate

    model_p = [ k0,  k1,  k2,  d,  m,  p,  k,  pp,  kk,  δ,  α, N];

    # println("params")
    # df_p = DataFrame();
    # df_p.variables = [:M,:R,:MR,:H4,:KDM5A,:H0,:PRC2,:H27,:KDM6A,:KMT];
    # df_p.values = model_p;
    # @show df_p
    # println("\n")

    ss = steady_states(Notch_model_cp, model_p)
    df = DataFrame(ss); df.vars = Notch_model_cp.syms
    @show df
    sort!(ss, by = x -> x[4])
    sb2 = stability_tianchi(ss,Notch_model_cp,model_p,1)
    # sb2 = stability(ss,Notch_model2,model_p)
    @show sb2
    println("\n")
end


