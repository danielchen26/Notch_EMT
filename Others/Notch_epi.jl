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
# using Interact, DiffEqBiological, DataFrames
# include(pwd()*"/Notch_Epi/functions.jl")
#
# Notch_model = @reaction_network begin
#
# # 	The competitve binding between NICD and MITF
# # 	model NICD as an external input
#
#   (N,k0),              R ↔ NR               # NICD binds RBPJ
#   (k1, k2),            M + R ↔ MR 		    # MITF binds RBPJ
#    # hill(M,v1,Kd1,2),   M --> M + R          # MITF form homodimer to activate RBPJ
#    # hillr(N,v2,Kd2,1),   N + M --> N         # NICD transcriptionally inhibit MITF
#
# # 	The epigenetic core regulation for the histone state of mir-222
#
# 	d1, 		H4  + KDM5A  --> H0  + KDM5A  # Demethylation of active mark
# 	m0,			H0  + PRC2   --> H27 + PRC2   # Methylation to get repressive mark
# 	d0, 		H27 + KDM6A  --> H0  + KDM6A  # Demethylation of repressive mark
# 	m1, 		H0  + KMT    --> H4  + KMT    # Methylation to get active mark
#
# # 	Epigenetic feeback
# 	k0,			MR --> MR + KDM5A    				# MITF-RBPJ recruit KDM5A
# 	H27, 		∅ --> PRC2           				# PRC2 is enhenced by H27 mark
# 	H4,	 		∅ --> KDM6A        					# KDM6A is enhenced by H4 mark
# 	H4, 		∅ --> KMT          					# KMT is enhenced by H4 mark
# 	H27, 		∅ --> KDM5A        				    # KDM5A is enhenced by H27 mark
# 	δ,			(PRC2, KDM5A, KDM6A, KMT) --> ∅     # Degradation of histone reader and writers
# # basal creation
# 	α1,          ∅ --> (KDM5A, KDM6A, KMT)
# 	α2,          ∅ --> PRC2
# end N k0 k1 k2 d1 d0 m1 m0 δ α1 α2 #v1 Kd1 v2 Kd2
# @add_constraints Notch_model  begin
#     R + NR + MR  = 10.0
# 	# H4 + H0 + H27 = 100.0
# end
#
# # default
# # model_p = [10.39, 2.8e-2, 2.0, 1.0,  9.13e-2,  8.39,  40.47,  8.39, 0.0,   1.0,  39.17,  1.0]
# # model_p = [N, k0, k1, k2, d1, d0, m1, m0, δ, v1, Kd1, v2, Kd2]
# model_p = [0.0, 10.0, 10.0, 10.0, 100.0, 100.0, 100.0, 100.0, 1.0, 1.0, 50.0, 1.0, 50.0, 1.0, 1.9]
# steady_states(Notch_model, model_p)


# @manipulate for
#     m1 = slider(0.01 : 0.01: 100.0, value=10.39),m2 = slider(0.01 : 0.01: 100.0, value=2.8e-2),
#     N = slider(0.00 : 0.01: 100.0, value= 2.0), r = slider(0.01 : 0.01: 100.0, value=1.0),
#     k0 = slider(0.01 : 0.01: 100.0, value= 0.09), k1 = slider(0.01 : 0.01: 100.0, value =8.93),
#     k2 = slider(0.01 : 0.01: 100.0, value= 40.47), k3 = slider(0.01 : 0.01: 100.0, value = 8.39),
#     k4 = slider(0.01 : 0.01: 100.0, value = 0.0),
#     d = slider(0.01 : 0.01: 100.0, value= 1.0),
#     p1 = slider(0.01 : 0.01: 100.0, value= 39.17), p2 = slider(0.01 : 0.01: 100.0, value= 1.0)
#
#     model_p = [   m1,      m2,   N,  r,     k0,      k1,     k2,    k3,  k4,      d,    p1,   p2]
#
#     println("params")
#     df_p = DataFrame()
#     df_p.variables = [:m1,:m2,:N,:r,:k0,:k1,:k2,:k3,:k4,:d,:p1,:p2];
#     df_p.values = model_p
#     @show df_p
#     println("\n")
#
#     ss = steady_states(NE_model3, model_p)
#     df = DataFrame(ss); df.vars = NE_model3.syms
#     @show df
#     sort!(ss, by = x -> x[3])
#     sb2 = stability_tianchi(ss,NE_model3,model_p,1)
#     @show sb2
#     println("\n")
# end


# #################################################################
# Use the following
#  #################################################################

# ## =================================================================
# using Interact, DiffEqBiological, DataFrames
# include(pwd()*"/Notch_Epi/functions.jl")
# Notch_model2 = @reaction_network begin
#   (k1, k2),            M + R ↔ MR 		    # MITF binds RBPJ
# 	d1, 		H4  + KDM5A  --> H0  + KDM5A  # Demethylation of active mark
# 	m0,			H0  + PRC2   --> H27 + PRC2   # Methylation to get repressive mark
# 	d0, 		H27 + KDM6A  --> H0  + KDM6A  # Demethylation of repressive mark
# 	m1, 		H0  + KMT    --> H4  + KMT    # Methylation to get active mark
#
# 	k0,			MR --> MR + KDM5A    				        # MITF-RBPJ recruit KDM5A
# 	p, 		    H27 --> H27 + PRC2           				# PRC2 is enhenced by H27 mark
# 	kk,	 		H4 --> H4 + KDM6A        					# KDM6A is enhenced by H4 mark
# 	pp, 		H4 --> H4 + KMT          					# KMT is enhenced by H4 mark
# 	k, 		    H27 --> H27 + KDM5A        				    # KDM5A is enhenced by H27 mark
# 	δ,			(PRC2, KDM5A, KDM6A, KMT) --> ∅             # Degradation of histone reader and writers
#
# 	α1,          ∅ --> (KDM5A, KDM6A, KMT)
# 	α2,          ∅ --> PRC2
# end k0 k1 k2 d1 d0 m1 m0 p k pp kk δ α1 α2
#
# @add_constraints Notch_model2  begin
#     M + MR  = 10.0
# 	R + MR  = 10.0
# 	H4 + H27 + H0 = 10.0
# end
#
# ##
# @show params(Notch_model2)
# model_p = [106.8, 137.62, 11.73, 0.913, 60.215, 27.57, 3.927, 223.89, 223.89, 27.57, 27.57, 1.0, 1.0, 1.98]
# ss = steady_states(Notch_model2, model_p)
# ##
# @manipulate for
#     k0 = slider(0.00 : 0.01: 100.0, value=106.8),
#     k1 = slider(0.00 : 0.01: 100.0, value= 137.62),
#     k2 = slider(0.00 : 0.01: 100.0, value= 11.73),
#     d1 = slider(0.00 : 0.01: 100.0, value= 0.913),
#     d0 = slider(0.00 : 0.01: 100.0, value = 60.215),
#     m1 = slider(0.00 : 0.01: 100.0, value= 27.57),
#     m0 = slider(0.00 : 0.01: 100.0, value= 3.927),
# 	p = slider(0.00 : 0.01: 100.0, value= 223.89),
# 	k = slider(0.00 : 0.01: 100.0, value=223.89),
# 	pp = slider(0.0 : 0.01: 100.0, value=27.57),
# 	kk = slider(0.0 : 0.01: 100.0, value = 27.57),
#  	δ = slider(0.0 : 0.01: 100.0, value = 1.00),
#  	α1 =  slider(0.0 : 0.01: 100.0, value= 1.00),
# 	α2 = slider(0.0 : 0.01: 100.0, value= 1.98)
#
#     model_p = [ k0,  k1,  k2,  d1,  d0,  m1,  m0,  p,  k,  pp,  kk,  δ,  α1,  α2];
#
#     # println("params")
#     # df_p = DataFrame();
#     # df_p.variables = [:M,:R,:MR,:H4,:KDM5A,:H0,:PRC2,:H27,:KDM6A,:KMT];
#     # df_p.values = model_p;
#     # @show df_p
#     # println("\n")
#
#     ss = steady_states(Notch_model2, model_p)
#     df = DataFrame(ss); df.vars = Notch_model2.syms
#     @show df
#     sort!(ss, by = x -> x[4])
#     sb2 = stability_tianchi(ss,Notch_model2,model_p,1)
#     # sb2 = stability(ss,Notch_model2,model_p)
#     @show sb2
#     println("\n")
# end



##
function make_cb(ts_in, index, value)
    ts = ts_in
    condition(u,t,integrator) = t in ts
    function affect!(integrator)
            if integrator.t == ts[1]
                integrator.p[index] = value
            elseif integrator.t == ts[2]
                integrator.p[index] = 0.0
            end
    end
    cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
    return ts, cb
end

"""
Randomly assign M,R,MR,H4,KDM5A,H0,PRC2,H27,KDM6A,KMT between 0 and 100, constrained by conservation law specified by H_c, M_c, R_c.
"""
function initialization_con(H_c, M_c, R_c)
	# inital conditions
	H_set = randfixsum(1, 3, H_c)
	M_con = M_c; R_con = R_c;
	MR = rand(0:minimum([M_c, R_c]),1)[1]
	NR = rand(0:minimum([M_c - MR, R_c - MR]),1)[1]
	M = M_con - MR
	R = R_con - MR - NR
	H4, H0, H27 = H_set[1][1], H_set[1][2], H_set[1][3]
	KDM5A, PRC2, KDM6A, KMT = rand(0:100,4)
	# u0 = [M,R,MR,H4,KDM5A,H0,PRC2,H27,KDM6A,KMT]
	u0 = [R,NR,M,MR,KDM5A,H4,H0,PRC2,H27,KDM6A,KMT]
end

# ## calculate the MR, H4, H27 equilibrium points
# result =[]
# u0_set = []
# for repeat = 1:10000
# 	# inital conditions
# 	u0 = initialization_con(50.0, 50.0, 50.0)
# 	push!(u0_set,u0)
# 	tspan = [0,1e2]
# 	prob = ODEProblem(Notch_rd_ODE,u0,tspan,model_p)
# 	sol = solve(prob)
# 	push!(result, sol[end][[3,4,8]])
# 	# push!(result, sol[end])
# end
# rr = [round.(i,digits = 4) for i in result]
# state_var = DataFrame(hcat(rr...)'); rename!(state_var,[:MR, :H4, :H27])
# u0_db = DataFrame(hcat(u0_set...)'); rename!(u0_db,[:u0_M, :u0_R, :u0_MR, :u0_H4, :u0_KDM5A, :u0_H0, :u0_PRC2, :u0_H27, :u0_KDM6A, :u0_KMT])
# db = hcat([state_var,u0_db]...)
# unique(rr)



## #################################
# The complete model
# I only want to  model the parameters of the relative methylation
####################################
using Interact, DiffEqBiological, DataFrames
include(pwd()*"/Notch_Epi/functions.jl")
# include(pwd()*"/functions.jl")
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

	# α,          ∅ --> (KDM5A, KDM6A, KMT, PRC2)
	α1,          ∅ --> (KDM6A, KMT, PRC2)
	α2,          ∅ --> (KDM5A )
# end k0 k1 k2 d m p k pp kk δ α N
end k0 k1 k2 d m p k pp kk δ α1 α2 N

@add_constraints Notch_model_cp  begin
    M + MR  = 50.0        # conservation of M
	R + MR + NR = 50.0     # conservation of R
	H4 + H27 + H0 = 50.0   # conservation of H
end

@show params(Notch_model_cp)
model_p = [1.0, 1.0, 1.0, 9.59, 1.05, 1.8, 3.77, 19.08, 19.08, 1.0, 1.98, 1.0, 0.0]
@show ss = steady_states(Notch_model_cp, model_p)


##
@manipulate for
    # k0 = slider(0.00 : 0.01: 100.0, value= 1.0,   label ="k0(MR induced KDM5A)"),     # MR --> KDM5A
    # k1 = slider(0.00 : 0.01: 100.0, value= 1.0, label ="k1(MR forming)"),   # M + R ↔ MR forward
    # k2 = slider(0.00 : 0.01: 100.0, value= 1.0,  label ="k2(MR diassembly)"),    # M + R ↔ MR backward
    # d = slider(0.00 : 0.0001: 10.0, value= 0.01,   label ="d(Demethylation ratio for H4)"),      # Demethylation of active mark
    # m = slider(0.00 : 0.01: 10.0, value= 2.52,   label ="m(Methylation ratio for H27)"),      # Methylation to get repressive mark
	# p = slider(0.00 : 0.01: 100.0, value= 19.89,  label ="p(Methy_H27)"),     # PRC2 is enhenced by H27 mark
	# k = slider(0.00 : 0.01: 100.0, value=19.89,   label ="k(Demethy_H4)"),      # KDM5A is enhenced by H27 mark
	# pp = slider(0.0 : 0.01: 100.0, value=2.57,    label ="pp(Methy_H4)"),      # KMT is enhenced by H4 mark
	# kk = slider(0.0 : 0.01: 100.0, value = 2.57,  label ="kk(Demthy_H27)"),    # KDM6A is enhenced by H4 mark
 	# δ = slider(0.0 : 0.01: 100.0, value = 1.00,    label ="δ(Epi degradation)"),      # Degradation of histone reader and writers
 	# # α =  slider(0.0 : 0.01: 100.0, value= 1.00,    label ="α(Epi production)"),      # Production rate
	# α1 =  slider(0.0 : 0.01: 100.0, value= 1.00,    label ="α1(Epi production)"),      # Production rate
	# α2 =  slider(0.0 : 0.01: 100.0, value= 1.00,    label ="α2(KDM5A production)"),      # Production rate
    # N =  slider(0.0 : 0.01: 100.0, value= 0.00,    label  ="N")     # Production rate
	k0 = slider(0.00 : 0.01: 100.0, value=1.0,   label ="k0(MR induced KDM5A)"),     # MR --> KDM5A
	k1 = slider(0.00 : 0.01: 100.0, value= 1.0, label ="k1(MR forming)"),   # M + R ↔ MR forward
	k2 = slider(0.00 : 0.01: 100.0, value= 1.0,  label ="k2(MR diassembly)"),    # M + R ↔ MR backward
	d = slider(0.00 : 0.0001: 10.02, value= 9.59,   label ="d(Demethylation ratio for H4)"),      # Demethylation of active mark
	m = slider(0.00 : 0.01: 10.0, value= 1.05,   label ="m(Methylation ratio for H27)"),      # Methylation to get repressive mark
	p = slider(0.00 : 0.01: 100.0, value= 1.8,  label ="p(Methy_H27)"),     # PRC2 is enhenced by H27 mark
	k = slider(0.00 : 0.01: 100.0, value= 3.77,   label ="k(Demethy_H4)"),      # KDM5A is enhenced by H27 mark
	pp = slider(0.0 : 0.01: 100.0, value=19.08,    label ="pp(Methy_H4)"),      # KMT is enhenced by H4 mark
	kk = slider(0.0 : 0.01: 100.0, value = 19.08,  label ="kk(Demthy_H27)"),    # KDM6A is enhenced by H4 mark
	δ = slider(0.0 : 0.01: 100.0, value = 1.98,    label ="δ(Epi degradation)"),      # Degradation of histone reader and writers
	α1 =  slider(0.0 : 0.01: 100.0, value= 1.00,    label ="α1(Epi production)"),      # Production rate
	α2 =  slider(0.0 : 0.01: 100.0, value= 1.00,    label ="α2(KDM5A production)"),      # Production rate
	N =  slider(0.0 : 0.01: 100.0, value= 0.00,    label  ="N")


    # model_p = [ k0,  k1,  k2,  d,  m,  p,  k,  pp,  kk,  δ,  α, N];
	model_p = [ k0,  k1,  k2,  d,  m,  p,  k,  pp,  kk,  δ,  α1,  α2, N]
    # println("params")
    # df_p = DataFrame();
    # df_p.variables = [:M,:R,:MR,:H4,:KDM5A,:H0,:PRC2,:H27,:KDM6A,:KMT];
    # df_p.values = model_p;
    # @show df_p
    # println("\n")
    ss = steady_states(Notch_model_cp, model_p)
    df = DataFrame(ss); df.vars = Notch_model_cp.syms
    @show df
    sort!(ss, by = x -> x[6])
    # sb2 = stability_tianchi(ss,Notch_model_cp,model_p,1)
    # sb2 = stability_tianchi(ss,Notch_model_cp,model_p,1)
    # @show sb2
    println("\n")
end

## get ub and lb for each parameters

for k0 = 0.00 : 1.0: 100.0,
    k1 = 0.00 : 0.01: 100.0,
    k2 = 0.00 : 0.01: 100.0,
    d = 0.00 : 0.0001: 0.02,
    m = 0.00 : 0.01: 10.0,
	p = 0.00 : 0.01: 100.0,
	k = 0.00 : 0.01: 100.0,
	pp = 0.0 : 0.01: 100.0,
	kk = 0.0 : 0.01: 100.0,
 	δ = 0.0 : 0.01: 100.0,
 	α =  0.0 : 0.01: 100.0,
    N =  0.0 : 0.01: 100.0
	model_p = [ k0,  k1,  k2,  d,  m,  p,  k,  pp,  kk,  δ,  α1,  α2, N]
	@show model_p
end




## calculate the MR, H4, H27 equilibrium points
function make_cb(ts_in, index, value)
    ts = ts_in
    condition(u,t,integrator) = t in ts
    function affect!(integrator)
            if integrator.t == ts[1]
                integrator.p[index] = value
            elseif integrator.t == ts[2]
                integrator.p[index] = 0.0
            end
    end
    cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
    return ts, cb
end

using DifferentialEquations

# modif this model with seperated activation rate of KDM5A
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

	α1,          ∅ --> (KDM6A, KMT, PRC2)
	α2,          ∅ --> (KDM5A)
end k0 k1 k2 d m p k pp kk δ α1 α2 N



@manipulate for
    # k0 = slider(0.00 : 0.01: 100.0, value=9.8,   label ="k0(MR induced KDM5A)"),     # MR --> KDM5A
    # k1 = slider(0.00 : 0.01: 100.0, value= 1.0, label ="k1(MR forming)"),   # M + R ↔ MR forward
    # k2 = slider(0.00 : 0.01: 100.0, value= 1.0,  label ="k2(MR diassembly)"),    # M + R ↔ MR backward
    # d = slider(0.00 : 0.0001: 0.02, value= 0.01,   label ="d(Demethylation ratio for H4)"),      # Demethylation of active mark
    # m = slider(0.00 : 0.01: 10.0, value= 2.52,   label ="m(Methylation ratio for H27)"),      # Methylation to get repressive mark
	# p = slider(0.00 : 0.01: 100.0, value= 19.89,  label ="p(Methy_H27)"),     # PRC2 is enhenced by H27 mark
	# k = slider(0.00 : 0.01: 100.0, value=19.89,   label ="k(Demethy_H4)"),      # KDM5A is enhenced by H27 mark
	# pp = slider(0.0 : 0.01: 100.0, value=2.57,    label ="pp(Methy_H4)"),      # KMT is enhenced by H4 mark
	# kk = slider(0.0 : 0.01: 100.0, value = 2.57,  label ="kk(Demthy_H27)"),    # KDM6A is enhenced by H4 mark
 	# δ = slider(0.0 : 0.01: 100.0, value = 1.00,    label ="δ(Epi degradation)"),      # Degradation of histone reader and writers
 	# α1 =  slider(0.0 : 0.01: 100.0, value= 1.00,    label ="α1(Epi production)"),      # Production rate
	# α2 =  slider(0.0 : 0.01: 100.0, value= 1.00,    label ="α2(KDM5A production)"),      # Production rate
    # N =  slider(0.0 : 0.01: 100.0, value= 0.00,    label  ="N")
	k0 = slider(0.00 : 0.01: 100.0, value=1.0,   label ="k0(MR induced KDM5A)"),     # MR --> KDM5A
	k1 = slider(0.00 : 0.01: 100.0, value= 1.0, label ="k1(MR forming)"),   # M + R ↔ MR forward
	k2 = slider(0.00 : 0.01: 100.0, value= 1.0,  label ="k2(MR diassembly)"),    # M + R ↔ MR backward
	d = slider(0.00 : 0.0001: 10.02, value= 9.59,   label ="d(Demethylation ratio for H4)"),      # Demethylation of active mark
	m = slider(0.00 : 0.01: 10.0, value= 1.05,   label ="m(Methylation ratio for H27)"),      # Methylation to get repressive mark
	p = slider(0.00 : 0.01: 100.0, value= 1.8,  label ="p(Methy_H27)"),     # PRC2 is enhenced by H27 mark
	k = slider(0.00 : 0.01: 100.0, value= 3.77,   label ="k(Demethy_H4)"),      # KDM5A is enhenced by H27 mark
	pp = slider(0.0 : 0.01: 100.0, value=19.08,    label ="pp(Methy_H4)"),      # KMT is enhenced by H4 mark
	kk = slider(0.0 : 0.01: 100.0, value = 19.08,  label ="kk(Demthy_H27)"),    # KDM6A is enhenced by H4 mark
	δ = slider(0.0 : 0.01: 100.0, value = 1.98,    label ="δ(Epi degradation)"),      # Degradation of histone reader and writers
	α1 =  slider(0.0 : 0.01: 100.0, value= 1.00,    label ="α1(Epi production)"),      # Production rate
	α2 =  slider(0.0 : 0.01: 100.0, value= 1.00,    label ="α2(KDM5A production)"),      # Production rate
	N =  slider(0.0 : 0.01: 100.0, value= 0.00,    label  ="N"),
	input =  slider(0.0 : 100000.0, value= 0.00,    label  ="Input")


#   test xor of before and after H27, existence means sum(row_select) > 0
	result =[]
	u0_set = []
	ti = 50.0; tf = 100.0; activation = input
	ts, cb  = make_cb([ti,tf], length(params(Notch_model_cp)), activation)
	sol_ti_set = []
	H27_set = []
	for repeat = 1:1000
		# inital conditions
		u0 = initialization_con(50.0, 50.0, 50.0)
		push!(u0_set,u0)
		tspan = [0,2e2]
		# model_p = [9.8, 1.0, 1.0, 0.01, 2.0, 20.0, 20.0, 2.57, 2.57, 1.0, 1.0, 0.0]
		model_p = [ k0,  k1,  k2,  d,  m,  p,  k,  pp,  kk,  δ,  α1,  α2, N];
		prob = ODEProblem(Notch_model_cp,u0,tspan,model_p)
		sol = solve(prob, callback=cb, tstops=ts)
		sol_ti = sol(ti)
		push!(sol_ti_set,sol_ti)
		push!(H27_set,sol(ti)[9])
		push!(result, sol[end][[4,5,6,9]])
		# push!(result, sol[end])
	end
	rr = [round.(i,digits = 2) for i in result]
	state_var = DataFrame(hcat(rr...)'); rename!(state_var,[:MR, :KDM5A, :H4, :H27])
	u0_db = DataFrame(hcat(u0_set...)'); rename!(u0_db,[:u0_R, :u0_NR, :u0_M, :u0_MR, :u0_KDM5A, :u0_H4, :u0_H0, :u0_PRC2, :u0_H27, :u0_KDM6A, :u0_KMT])
	db = hcat([state_var,u0_db]...)
	db.before_activation_H27 = H27_set .> 25.0
	db.after_activation_H27 = db.H27 .> 25.0
	db

	row_select = (db[:,end-1] .⊻  db[:,end]) .& db.before_activation_H27 # select the row
	# @show row_select
	@show db[row_select,:]
	# @show sum(row_select)

println("\n")
end


# state_var_count = combine(groupby(state_var, [:MR, :KDM5A, :H4, :H27]), nrow => :count)
# @show state_var_count




## ---------
using Plots;gr()
ts, cb  = make_cb([50.0,100.0], 12, 1000.0)
# u0 = initialization_con(50.0, 500.0, 50.0)
# u0 = [17.0, 6.0, 473.0,  27.0,  21.0,  13.7,  23.8,  37.0,  12.4,  12.0,  37.0] # repressive mark
u0 = [11.0    ,3.0     ,464.0   ,36.0    ,48.0     ,13.6876 ,10.678  ,38.0    ,25.6344 ,57.0     ,13.0    ]
tspan = (0.0,150)
model_p = [9.8, 1.0, 1.0, 0.001, 2.0, 20.0, 20.0, 2.57, 2.57, 1.0, 1.0, 0.0]
prob = ODEProblem(Notch_model_cp, u0, tspan, model_p)
sol = solve(prob, callback=cb, tstops=ts)
plot(sol, vars = [ :H4, :H27], lw  = 1.5)#:MR, :KDM5A,

u0_set = db[:,[:u0_R, :u0_NR, :u0_M, :u0_MR, :u0_KDM5A, :u0_H4, :u0_H0, :u0_PRC2, :u0_H27, :u0_KDM6A, :u0_KMT]]
for i = 1:10
	u0 = Array(eachrow(u0_set)[i])
	@show u0
	prob = ODEProblem(Notch_model_cp, u0, tspan, model_p)
	sol = solve(prob, callback=cb, tstops=ts)
	plt = plot(sol, vars = [ :H4, :H27], lw  = 1.5, title = "$i")#:MR, :KDM5A,
	display(plt)
end









## test the core model
core_model = @reaction_network begin
# 	The epigenetic core regulation for the histone state of mir-222
	d, 		H4  + KDM5A  --> H0  + KDM5A        # Demethylation of active mark
	m,		H0  + PRC2   --> H27 + PRC2         # Methylation to get repressive mark
	1, 		H27 + KDM6A  --> H0  + KDM6A        # Demethylation of repressive mark
	1, 		H0  + KMT    --> H4  + KMT          # Methylation to get active mark
# 	Epigenetic feeback
	p, 		    H27 --> H27 + PRC2    				# PRC2 is enhenced by H27 mark
	kk,	 		H4  --> H4  + KDM6A   				# KDM6A is enhenced by H4 mark
	pp, 		H4  --> H4  + KMT     				# KMT is enhenced by H4 mark
	k, 		    H27 --> H27 + KDM5A   				# KDM5A is enhenced by H27 mark
	δ,			(PRC2, KDM5A, KDM6A, KMT) --> ∅     # Degradation of histone reader and writers
# 	Basal production rates for epigenetic writers
	α1,          ∅ --> (KDM6A, KMT, PRC2)
	α2,			 ∅ --> (KDM5A)

end d m p k pp kk δ α1 α2 k0

@add_constraints core_model  begin
	H4 + H27 + H0 = 50.0   # conservation of H
end

@show params(core_model)
model_p = [ 1.0, 1.0, 1.0, 9.59, 1.05, 1.8, 3.77, 19.08, 19.08, 10.0]
@show ss = steady_states(core_model, model_p)


##

@manipulate for
	d = slider(0.00 : 0.0001: 10.02, value= 9.59,   label ="d(Demethylation ratio for H4)"),      # Demethylation of active mark
	m = slider(0.00 : 0.01: 10.0, value= 1.05,   label ="m(Methylation ratio for H27)"),      # Methylation to get repressive mark
	p = slider(0.00 : 0.01: 100.0, value= 1.8,  label ="p(Methy_H27)"),     # PRC2 is enhenced by H27 mark
	k = slider(0.00 : 0.01: 100.0, value= 3.77,   label ="k(Demethy_H4)"),      # KDM5A is enhenced by H27 mark
	pp = slider(0.0 : 0.01: 100.0, value=19.08,    label ="pp(Methy_H4)"),      # KMT is enhenced by H4 mark
	kk = slider(0.0 : 0.01: 100.0, value = 19.08,  label ="kk(Demthy_H27)"),    # KDM6A is enhenced by H4 mark
	δ = slider(0.0 : 0.01: 100.0, value = 1.98,    label ="δ(Epi degradation)"),      # Degradation of histone reader and writers
	α1 =  slider(0.0 : 0.01: 100.0, value= 1.00,    label ="α1(Epi production)"),      # Production rate
	α2 =  slider(0.0 : 0.01: 100.0, value= 1.00,    label ="α2(KDM5A production)")     # Production rate

	model_p = [d, m, p, k, pp, kk, δ, α1, α2]
	ss = steady_states(core_model, model_p)
	rr = [round.(i,digits = 1) for i in ss]
	df = DataFrame(rr); df.vars = core_model.syms
	@show df
	sort!(ss, by = x -> x[6])
	println("\n")
end







##    analytical
import SymPy

SymPy.@vars R NR M MR k1 k2 N
SymPy.@vars RT MT

@show M_ss =SymPy.solve(k2*(M - MT) + M*k1*(M - MT + RT),M)
