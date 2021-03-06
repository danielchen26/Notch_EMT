##
using DifferentialEquations, DiffEqBiological
using DataFrames, LinearAlgebra, ProgressMeter
# using Plots


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
    M + MR  = 50.0        # conservation of M
	R + MR + NR = 50.0     # conservation of R
	H4 + H27 + H0 = 50.0   # conservation of H
end

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

# tianchi_stability_local(ss, Notch_model_cp, model_p)

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

# # One instance of example for model_p
model_p = [6.0, 1.0, 1.0, 0.1, 0.4, 5.0, 4.0, 4.0, 4.0, 1.0, 1.0, 0.0]
steady_states(Notch_model_cp, model_p)

# @showprogress for  k0 = 1.0:4.0: 10.0,  k1 = 1.0,  k2 = 1.0,
#      d = 0.01:0.2:1.0,  m = 0.01:0.2:1.0,  p = 1.0:5:20.0,
#      k = 0.0:5: 20.0,  pp = 1.0:5:20.0,  kk = 1.0:5:20.0,
#      δ = 1.0,  α1 = 1.0, N = 0.0
@showprogress for  k0 = 6.0,  k1 = 1.0,  k2 = 1.0,
      d = 0.1,  m = 0.4,  p = 5.0,
      k = 0.0:4: 20.0,  pp = 0.0:4: 20.0,  kk = 0.0:4: 20.0,
      δ = 1.0,  α1 = 1.0, N = 0.0 # this is a short working example

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
		# @show  ss_H27_high
		# start from the ss_H27_high state
		prob = ODEProblem(Notch_model_cp,ss_H27_high,(0.0, t_end),model_p)
		sol = solve(prob, callback=cb, tstops=ts)
		# sol = solve(prob)
		if  sol(ti - 5.0)[9] > sol(ti - 5.0)[6] && sol(t_end)[6] > sol(t_end)[9]
			# plt = plot(sol, vars = [ :MR, :KDM5A, :H4, :H27], lw  = 1.5, title = "$model_p")
			# display(plt)
			# record the paramters to a DataFrame
			push!(parameters_set, model_p)
			push!(H4_H27_states, [sol(ti - 5.0)[6], sol(ti - 5.0)[9], sol(t_end)[6], sol(t_end)[9]])
			push!(initial_H27_high_states, ss_H27_high)
		end
	end
	# if tt >= 100 # plot 10 examples
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

using CSV
CSV.write("Notch_params_H4H27.csv", database_H4H27)
CSV.write("Notch_params_complete.csv", database_complete)
