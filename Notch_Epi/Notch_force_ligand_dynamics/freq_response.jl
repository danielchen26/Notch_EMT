using DifferentialEquations, Plots
using ModelingToolkit, Catalyst
# using
using Interact,ProgressMeter
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
##
model = @reaction_network begin
    # (A*(1+sign(cos(w*t))),1.0),     R ↔ NR               					# NICD binds RBPJ
    (A*(1+cos(w*t)),1.0),     R ↔ NR               					# NICD binds RBPJ
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
	α1,         ∅ --> (KDM6A, KMT, PRC2, KDM5A)
end k0 k1 k2 d m p k pp kk δ α1 w A # put A at last as the control 12th variable
@show speciesmap(model)
@show paramsmap(model)


u0 = [6.0, 0.0, 6.0, 40.0, 500.0, 0.6876,  0.678,  500.0,  50.6344, 1.0, 2.0]
tspan = [0.0, 150.0]

# # with w =0, when N is turned on, it is a constant rate
# p = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, 0.0, 0.0]

# with w = 1,when N is turned on, it is a oscillating rate
p = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, 0.5, 0.0]

ts, cb  = make_cb([50.0,100.0], 13, 80.0)
prob = ODEProblem(model, u0, tspan, p)
@time sol = solve(prob, Rosenbrock23(), callback=cb, tstops=ts)
plot(sol, vars = [2,5,6,9,10], lw  = 1.5)


## frequency response
