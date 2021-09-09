using Interact 
using Catalyst, DataFrames
using OrdinaryDiffEq, Plots;plotly()
using Latexify, HomotopyContinuation



Notch_model_cp = @reaction_network begin
    (N,1.0),              R ↔ NR               					# NICD binds RBPJ
    (k1, k2),             M + R ↔ MR 		 			        # MITF binds RBPJ
     k0,		 		 MR --> MR + KDM5A    				        	# MITF-RBPJ recruit KDM5A
	d, 	      	 		H4  + KDM5A  --> H0  + KDM5A      				 # Demethylation of active mark
	m,					H0  + PRC2   --> H27 + PRC2   					# Methylation to get repressive mark
	1.0, 				H27 + KDM6A  --> H0  + KDM6A 				 # Demethylation of repressive mark
	1.0, 				H0  + KMT    --> H4  + KMT   					 # Methylation to get active mark
	p, 			H27 --> H27 + PRC2           					# PRC2 is enhenced by H27 mark
	kk,	 		H4 --> H4 + KDM6A        						# KDM6A is enhenced by H4 mark
	pp, 		H4 --> H4 + KMT          						# KMT is enhenced by H4 mark
	k, 			H27 --> H27 + KDM5A        				        # KDM5A is enhenced by H27 mark
	δ,			(PRC2, KDM5A, KDM6A, KMT) --> ∅            		        # Degradation of histone reader and writers
	α1,          ∅ --> (KDM6A, KMT, PRC2, KDM5A)
end k0 k1 k2 d m p k pp kk δ α1 N # put N at last as the control 12th variable



# the ode model of the above crn, 
@var R NR M MR KDM5A H4 H0 PRC2 H27 KDM6A KMT
function steady_states(param)
    k0, k1, k2, d, m, p, k, pp, kk, δ, α1, A = param
    @show k0, k1, k2, d, m, p, k, pp, kk, δ, α1, A
    f_R = -A*R + NR - k1*M*R + k2*MR
    f_NR = A*R - NR
    f_M = -k1*M*R + k2*MR 
    f_MR = k1*M*R - k2*MR 
    f_H4 = -d*H4*KDM5A + H0*KMT
    f_KDM5A = k0*MR + k*H27 -δ*KDM5A + α1
    f_H0 = d*H4*KDM5A - m*H0*PRC2 + H27*KDM6A - H0*KMT
    f_PRC2 = p*H27 - δ*PRC2 + α1
    f_H27 = m*H0*PRC2 - H27*KDM6A 
    f_KDM6A = kk*H4 - δ*KDM6A + α1
    f_KMT = pp*H4 - δ*KMT + α1 
    #conservation law
    Con_M = M + MR - 50.0   # conservation of M
	Con_R = R + MR + NR - 50.0     # conservation of R
	Con_H = H4 + H27 + H0 - 50.0
    F = System([f_R, f_NR, f_M, f_MR, f_H4, f_KDM5A, f_H0, f_PRC2, f_H27, f_KDM6A, f_KMT, Con_M, Con_R, Con_H])
    result = HomotopyContinuation.solve(F)
    real_solution = real_solutions(result)
    idx = [all(>=(0), i) for i in real_solution]
    positive_sol = real_solution[idx]
end

model_p = [6.0, 1.0, 1.0, 0.1, 0.4, 5.0, 4.0, 4.0, 4.0, 1.0, 1.0, 0.0]
k0, k1, k2, d, m, p, k, pp, kk, δ, α1, A = model_p
positive_sol = steady_states(model_p)