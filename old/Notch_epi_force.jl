using Interact, Latexify, DataFrames, Catalyst
using DifferentialEquations, Plots
using ModelingToolkit


@variables R NR M MR KDM5A H4 H0 PRC2 H27 KDM6A KMT

##
# @parameters w A
# @variables t x(t) f(t)
# D = Differential(t)

# @named test_model = ODESystem([f ~ A*cos(w*t), D(x) ~ f])

# # with out callback
# prob = ODEProblem(structural_simplify(test_model), [x => 1.0], (0.0,100.0),[w => 1.0, A => 1.0])
# sol = solve(prob)
# plot(sol, vars=[x,f])




# with callback controling parameter A
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

ts, cb  = make_cb([50.0,80.0], 1, 2.0)
prob = ODEProblem(structural_simplify(test_model), [x => 1.0], (0.0,100.0),[w => 1.0, A => 0.0])
sol = solve(prob, callback=cb, tstops=ts)
plot(sol, vars=[x,f])







## The notch epi model with oscillatory Notch ligand dynamics

function make_cb_notch(ts_in, index, value)
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


Notch_model_cp = @reaction_network begin
    (A*(1+cos(w*t)),1.0),              R ↔ NR               					# NICD binds RBPJ
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
end k0 k1 k2 d m p k pp kk δ α1 w A # put N at last as the control 12th variable

@show paramsmap(Notch_model_cp)
@show speciesmap(Notch_model_cp)
odesys = convert(ODESystem, Notch_model_cp)

ts, cb  = make_cb_notch([50.0,100.0], 13, 1000.0)

u0 = [6.0    ,0.0     ,6.0   ,40.0    ,500.0     ,0.6876 ,0.678  ,500.0    ,50.6344 ,1.0     ,2.0    ]
tspan = (0.0,150)
model_p = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, 1.0, 0.0]

u0map  = Pair.(species(Notch_model_cp), u0)
pmap   = Pair.(params(Notch_model_cp), model_p)
prob2 = ODEProblem(odesys, u0map, tspan, pmap)

# prob = ODEProblem(Notch_model_cp, u0, tspan, model_p)
sol = solve(prob2, callback=cb, tstops=ts)
dplt = plot(sol, vars = [4,5,6,9], lw  = 1.5)


