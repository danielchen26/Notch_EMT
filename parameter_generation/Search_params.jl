using SymPy
@vars R NR M MR H4 KDM5A H0 PRC2 H27 KDM6A KMT
@vars k0 k1 k2 d m p k pp kk δ α1 N

eqs = [ NR - N*R + k2*MR - M*R*k1
        -NR + N*R
        k2*MR - M*R*k1
        -k2*MR + M*R*k1
        H0*KMT - d*H4*KDM5A
        α1 + k*H27 + k0*MR - δ*KDM5A
        -H0*KMT + H27*KDM6A + d*H4*KDM5A - m*H0*PRC2
        α1 + p*H27 - δ*PRC2
        -H27*KDM6A + m*H0*PRC2
        α1 + kk*H4 - δ*KDM6A
        α1 + pp*H4 - δ*KMT ]

vars = [R, NR, M, MR, H4, KDM5A, H0, PRC2, H27, KDM6A, KMT]
# fps = SymPy.solve(eqs, vars)
J = eqs.jacobian(vars)



## Define model
using OrdinaryDiffEq, DiffEqBiological
using DataFrames, LinearAlgebra
using Plots
using ProgressMeter
# include(pwd()*"/Notch_Epi/functions.jl")


Notch_model_cp = @reaction_network begin
    (N,1.0),              R ↔ NR               					# NICD binds RBPJ
    (k1, k2),             M + R ↔ MR 		 			        # MITF binds RBPJ
     k0,		  MR --> MR + KDM5A    				        # MITF-RBPJ recruit KDM5A
	d, 	       	H4  + KDM5A  --> H0  + KDM5A  # Demethylation of active mark
	m,		H0  + PRC2   --> H27 + PRC2   # Methylation to get repressive mark
	1.0, 		H27 + KDM6A  --> H0  + KDM6A  # Demethylation of repressive mark
	1.0, 		H0  + KMT    --> H4  + KMT    # Methylation to get active mark
	p, 		H27 --> H27 + PRC2           					# PRC2 is enhenced by H27 mark
	kk,	 	H4 --> H4 + KDM6A        					# KDM6A is enhenced by H4 mark
	pp, 		H4 --> H4 + KMT          					# KMT is enhenced by H4 mark
	k, 		H27 --> H27 + KDM5A        				        # KDM5A is enhenced by H27 mark
	δ,		(PRC2, KDM5A, KDM6A, KMT) --> ∅            		        # Degradation of histone reader and writers
	α1,          ∅ --> (KDM6A, KMT, PRC2, KDM5A)
end k0 k1 k2 d m p k pp kk δ α1 N # put N at last as the control 12th variable

@add_constraints Notch_model_cp  begin
    M + MR  = 50.0                      # conservation of M
	R + MR + NR = 50.0                  # conservation of R
	H4 + H27 + H0 = 50.0                # conservation of H
end

@show paramsmap(Notch_model_cp)
@show speciesmap(Notch_model_cp)
model_p = [1.0, 1.0, 1.0, 9.59, 1.05, 1.8, 3.77, 19.08, 19.08, 1.0, 1.98, 0.0]
@show ss = steady_states(Notch_model_cp, model_p)

df = DataFrame(ss)


## Substitute fixed point values
idx = 3
J_fp = J.subs([(R,ss[idx][1]),(NR,ss[idx][2]),(M,ss[idx][3]),(MR,ss[idx][4]),(KDM5A,ss[idx][5]),(H4,ss[idx][6]),(H0,ss[idx][7]),(PRC2,ss[idx][8]),(H27,ss[idx][9]),(KDM6A,ss[idx][10]),(KMT,ss[idx][10])])
J_eval = J_fp.subs([(k0,model_p[1]),(k1,model_p[2]),(k2,model_p[3]),(d,model_p[4]),(m,model_p[5]),(p,model_p[6]),(k,model_p[7]),(pp,model_p[8]),(kk,model_p[9]),(δ,model_p[10]),(α1,model_p[11]),(N,model_p[12])])

@time J_vales =map(x->convert(Float64,x),J_eval)
eigs = eigvals(J_vales)


## Find the stability of the fixed point
is_stable = false
if maximum(real(eigs)) < 0
    is_stable = true
end
is_stable



## Calculate the steady states solution and the st

"""
	Generate a dataframe of the steady state solution for the model
"""
function SS_sol(model, parameters)
	@show ss = steady_states(model, parameters)
	ss = [round.(i,digits = 2) for i in ss]
	df = DataFrame(ss)
end


df_ss = SS_sol(Notch_model_cp,model_p)



"""
	Calculate the Jacobian of the model and return the eigenvalues to determine the stability of the steady state solution
"""
function stability_local(eqs, vars, SS_sol, model_p)
	J = eqs.jacobian(vars)
	@show SS_sol

	# calculate the eigenvalues of the jacobian matrix for each steady state
	eigs_set = []
	for idx = 1 : ncol(SS_sol)
		J_fp = J.subs([(R,SS_sol[!, idx][1]),(NR,SS_sol[!, idx][2]),(M,SS_sol[!, idx][3]),(MR,SS_sol[!, idx][4]),(KDM5A,SS_sol[!, idx][5]),(H4,SS_sol[!, idx][6]),(H0,SS_sol[!, idx][7]),(PRC2,SS_sol[!, idx][8]),(H27,SS_sol[!, idx][9]),(KDM6A,SS_sol[!, idx][10]),(KMT,SS_sol[!, idx][10])])
		J_eval = J_fp.subs([(k0,model_p[1]),(k1,model_p[2]),(k2,model_p[3]),(d,model_p[4]),(m,model_p[5]),(p,model_p[6]),(k,model_p[7]),(pp,model_p[8]),(kk,model_p[9]),(δ,model_p[10]),(α1,model_p[11]),(N,model_p[12])])
		J_values =map(x->convert(Float64,x),J_eval)
		eigs = eigvals(J_values)
		@show eigs
		push!(eigs_set, eigs)
	end

	# calculate the local stability of each steady state
	is_stable_set = []
	for idx = 1 : ncol(SS_sol)
		if maximum(real(eigs_set[idx])) < 1e-7
		    is_stable = true
		    push!(is_stable_set, is_stable)
	    	else maximum(real(eigs_set[idx])) > 1e-7
		    is_stable = false
		    push!(is_stable_set, is_stable)
		end
	end
	return eigs_set, is_stable_set
end


eigs_set, is_stable_set = stability_local(eqs, vars, df_ss, model_p)




##
# Define local stability function
function tianchi_stability_local(ss, model, p)

	function JE_stability(solution::Vector{Float64}, rn::DiffEqBase.AbstractReactionNetwork, p::Vector{Float64}, t=0.::Float64)
	    jac = zeros(length(rn.syms),length(rn.syms))
	    rn.jac(jac,solution,p,t)
	    return (jac,eigen(jac).values)
	end
	Eigen_spectrum = [JE_stability(i, model, p)[2] for i in ss]
	tianchi_ss =  [maximum(real(i)) < 1e-10 for i in Eigen_spectrum]
end

tianchi_stability_local(ss, Notch_model_cp, model_p)

# test switch with callbacks starting from neat H27 state
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
ti = 50.0; tf = 100.0; t_end = 150.0; activation = 1000.0;
ts, cb  = make_cb([ti,tf], length(params(Notch_model_cp)), activation)


tt = 0
parameters_set = []
H4_H27_states = []
initial_H27_high_states = []
@showprogress for  k0 = 1.0:4.0: 30.0,  k1 = 1.0,  k2 = 1.0,
     d = 0.01:0.01:1.0,  m = 0.01:0.1:1.0,  p = 1.0:2:20.0,
     k = 0.0:20.0,  pp = 1.0:2:20.0,  kk = 1.0:2:20.0,
     δ = 1.0,  α1 = 1.0, N = 0.0
 # @showprogress for  k0 = 1.0:4.0: 30.0,  k1 = 1.0,  k2 = 1.0,
 #      d = 0.0:0.01:1.0,  m = 0.0:0.1:1.0,  p = 5.0,
 #      k = 1.0,  pp = 2.57,  kk = 8.0,
 #      δ = 1.0,  α1 = 1.0, N = 0.0 # this is a short working example

	model_p = [k0, k1, k2, d, m, p, k, pp, kk, δ, α1, N]
	@show model_p
	ss = steady_states(Notch_model_cp, model_p)
	sort!(ss, by = x -> x[9]) # assending H27
	# @show DataFrame(ss)
	stability = tianchi_stability_local(ss, Notch_model_cp, model_p)
	# @show stability
	# println("\n")

	if sum(stability) >= 2 && ss[end][9] > ss[end][6]# at least 2 stable fixed points
		ss_H27_high = ss[end]
		@show  ss_H27_high
		# start from the ss_H27_high state
		prob = ODEProblem(Notch_model_cp,ss_H27_high,(0.0, t_end),model_p)
		sol = solve(prob, callback=cb, tstops=ts)
		# sol = solve(prob)
		if  sol(ti - 5.0)[9] > sol(ti - 5.0)[6] && sol(t_end)[6] > sol(t_end)[9]
			plt = plot(sol, vars = [ :MR, :KDM5A, :H4, :H27], lw  = 1.5, title = "$model_p")
			display(plt)
			# record the paramters to a DataFrame
			push!(parameters_set, model_p)
			push!(H4_H27_states, [sol(ti - 5.0)[6], sol(ti - 5.0)[9], sol(t_end)[6], sol(t_end)[9]])
			push!(initial_H27_high_states, ss_H27_high)
		end
	end
	# if tt >= 3000 # plot 10 examples
	#       break
	# end
	# tt+=1
end


# @show [i[9] > i[6] for i in ss]

##
function db_gen(parameters_set, H4_H27_states, initial_H27_high_states)
	params_db = DataFrame(hcat(parameters_set...)')
	rename!(params_db,  [:k0, :k1, :k2, :d, :m, :p, :k, :pp, :kk, :δ, :α1, :N])

	# showall(params_db)

	H4_H27_db = DataFrame(hcat(H4_H27_states...)')
	rename!(H4_H27_db, [:H4_bf_ti, :H27_bf_ti, :H4_end, :H27_end])
	H4_H27_db

	initial_state_db = DataFrame(hcat(initial_H27_high_states...)')
	rename!(initial_state_db, [:R_init ,:NR_init ,:M_init ,:MR_init ,:KDM5A_init ,:H4_init ,:H0_init ,:PRC2_init ,:H27_init ,:KDM6A_init ,:KMT_init])

	database_H4H27 = hcat(params_db, H4_H27_db)
	database_complete = hcat(params_db,initial_state_db)
	return database_H4H27, database_complete
end

database_H4H27, database_complete = db_gen(parameters_set, H4_H27_states, initial_H27_high_states)

##
using CSV
CSV.write("Notch_params_H4H27_finer.csv", database_H4H27)
CSV.write("Notch_params_complete_finer.csv", database_complete)
