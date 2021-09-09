### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ dfb1f4ba-7ded-4b87-ac7d-4974349ed12d
using Catalyst, DifferentialEquations, Plots, PlutoUI, Images, Latexify

# ╔═╡ 8ea69aeb-13ef-4570-9150-4e23219fe8d4
using DataFrames

# ╔═╡ 3d428d36-cdf9-4f2b-99ee-14e7cbc1f39b
using DiffEqSensitivity, Statistics

# ╔═╡ 347867d4-726a-4e83-977a-fd5e77512422
md"""
# Notch Epigenetics
"""

# ╔═╡ d72bc2d4-b13f-450e-a24d-a056d7b48199
# TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ 42fb32c8-10a6-4e75-940b-07193d4e549b
import DiffEqBiological:steady_states

# ╔═╡ b5302649-d9be-47e3-9a11-d39370bb678a
plotly()

# ╔═╡ 968fa436-802a-4ff8-bc9a-f225ac439e9d
html"""<style>
main {
    max-width: 970px;
}
"""

# ╔═╡ 9e4cd586-71c0-4e7e-b958-1efe64f228b8


# ╔═╡ 2cc71eb3-df4f-4cf9-9682-1938959d58a2


# ╔═╡ d8346212-0526-4423-be21-74cde773322d
md"""
## Defining  the reduced model
The reduced model tend to directly modeling the dynamics of the epigenetic regulation which involves __active histone mark__ `H3K4me3`, __repressive histone mark__`H3K27me3`, and the epigenetic regulators `KDM5A`, `KDM6A`, `PRC2`, `KMT`.
"""

# ╔═╡ 528975db-a04c-4e2c-8aee-be2083ea6056
model_graph = "/Users/chentianchi/Desktop/Projects/Notch_EMT/Pluto_nb/pic/notch_EMT_graph.png"

# ╔═╡ 825f90df-6828-440e-b7b5-07106c1ae788
load(model_graph)

# ╔═╡ 36a50cf4-2b06-4e4b-9e66-469753e4f3b4
md"""
### Model variables
- N : `NICD`, **the introcellar domain of Notch receptor**.
- R : `RBPJ`, a DNA binding protein that act as either repressor or activator depands on the which cofactor it binds to.
- M : `MITF`, a protein that binds to `RBPJ`
- MR: The `MITF-RBPJ` complex, which will recruit the histone demethylate `KDM5A` to the promoter region of the mir-222 to get rid of the active histone mark `H3K4me3`.
- NR: `NICD-RBPJ` complex
- H4, the __active__ promoter state of mir-222 with histone mark `H3K4me3`
- H0, the __inactive__ promoter state of mir-222 without histone mark.
- H27, the __repressive__ promoter state of mir-222 with histone mark `H3k27me3`.
"""

# ╔═╡ 8732d482-987a-40a8-9da7-8ec9523931bb
md"""
### Chemical reaction network (CRN)
"""

# ╔═╡ e7d396f4-ecf6-4d33-80ff-2fac41990906
md"""
### From CRN to ODE
"""

# ╔═╡ 35028063-80b0-4ab8-96b4-bc84af7389b9
md"""
Define the function `initialization_con` to sample initial conditions for u0 with a fixed set of conservation laws. Namely:
- M + MR  = constant
- R + MR + NR  = constant
- H4 + H27 + H0 = constant
"""

# ╔═╡ e86a64f0-1997-4149-b87a-e55b909cc153
md"""
Define `initialization_con` and `randfixsum`
"""

# ╔═╡ 5bb6836f-698f-4276-a4bd-84c3260febfb


# ╔═╡ a44a95ac-7cb7-4b4b-b64d-7dd43ec09629
begin
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
		KDM5A, PRC2, KDM6A, KMT = rand(0:30,4)
		# u0 = [M,R,MR,H4,KDM5A,H0,PRC2,H27,KDM6A,KMT]
		u0 = [R,NR,M,MR,H4,KDM5A,H0,PRC2,H27,KDM6A,KMT]
	end


	"""
	Random numbers with fixed sum function
	"""
	function randfixsum(n, dim, tot)
	    sol = []
	    for i in 1: n
	        v = rand(1,dim)
	        nv = tot*v./sum(v)
	        push!(sol, nv)
	    end
	    return sol
	end
end

# ╔═╡ b5078a8f-8a6c-4257-ba4a-0f3ca0fac4c0
# begin
	
# 	u0 = initialization_con(50.0, 500.0, 50.0)
# 	tspan = [0,1e2]
# 	model_p = [k0, k1, k2, d, m, p, k, pp, kk, δ, α, N]
# 	prob = ODEProblem(Notch_model, u0, tspan, model_p)
# 	sol  = solve(prob);
# 	plot(sol, vars = [:MR, :KDM5A, :H4, :H27])
	
# end

# ╔═╡ 1eeae429-3c92-4bbf-9856-c7db46ccb95b
md"""
## Intercept the model by controling N
"""

# ╔═╡ 3fef48e3-a6f9-4c64-bdea-41a750c4e6bf
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

# ╔═╡ d47a550f-473e-47a7-bd89-e273f8f0854f


# ╔═╡ a360a666-653d-4aed-b08d-88bdb955773a
begin
	md"""
	### Time and activation strength
	tmax = $(@bind tmax Slider(0.0:1000.0, show_value=true, default = 150)),\
	ti = $(@bind ti Slider(0.0:1000.0, show_value=true, default = 50)),\
	tf = $(@bind tf Slider(0.0:1000.0, show_value=true, default = 80)),\
	activation = $(@bind activation Slider(0.0:10000.0, show_value=true, default = 1000))
	"""
end

# ╔═╡ 2b189db9-708a-49a1-8ef1-cc4a360fa8db
# begin
# 	md"""
# 	### Model parameters:
# 	k0 = $(@bind k0 Slider(0.0:50.0, show_value=true, default = 9.8)),
# 	k1 = $(@bind k1 Slider(0.0:50.0, show_value=true, default = 1.)),
# 	k2 = $(@bind k2 Slider(0.0:50.0, show_value=true, default = 1.)),\
# 	d = $(@bind d Slider(0.0: 0.01 : 10.01, show_value=true, default = 9.59)),
# 	m = $(@bind m Slider(0.0: 0.01 : 10.0, show_value=true, default = 1.05)),
# 	p = $(@bind p Slider(0.0:50.0, show_value=true, default = 1.8)),\
# 	k = $(@bind k Slider(0.0:50.0, show_value=true, default = 3.77)),
# 	pp = $(@bind pp Slider(0.0:50.0, show_value=true, default = 19.08)),
# 	kk = $(@bind kk Slider(0.0:50.0, show_value=true, default = 19.08)),\
# 	δ = $(@bind δ Slider(0.0:50.0, show_value=true, default = 1.98)),
# 	α1 = $(@bind α1 Slider(0.0:0.01:10.0, show_value=true, default = 1.0)),
# 	α2 = $(@bind α2 Slider(0.0:0.01:10.0, show_value=true, default = 1.0)),
# 	N = $(@bind N Slider(0.0:50.0, show_value=true, default = 0.0)),
# 	"""
# end
begin
	md"""
	### Model parameters(new):
	k0 = $(@bind k0 Slider(0.0:100.0, show_value=true, default = 100.)),
	k1 = $(@bind k1 Slider(0.0:50.0, show_value=true, default = 1.)),
	k2 = $(@bind k2 Slider(0.0:50.0, show_value=true, default = 1.)),\
	d = $(@bind d Slider(0.0: 0.01 : 10.01, show_value=true, default = 0.01)),
	m = $(@bind m Slider(0.0: 0.01 : 10.0, show_value=true, default = .51)),
	p = $(@bind p Slider(0.0:50.0, show_value=true, default = 5.0)),\
	k = $(@bind k Slider(0.0:50.0, show_value=true, default = 1.0)),
	pp = $(@bind pp Slider(0.0:50.0, show_value=true, default = 2.57)),
	kk = $(@bind kk Slider(0.0:50.0, show_value=true, default = 8.0)),\
	δ = $(@bind δ Slider(0.0:50.0, show_value=true, default = 5.0)),\
	α1 = $(@bind α1 Slider(0.0:0.01:10.0, show_value=true, default = 0.1)),
	α2 = $(@bind α2 Slider(0.0:0.01:10.0, show_value=true, default = 1.0)),
	α3 = $(@bind α3 Slider(0.0:0.01:10.0, show_value=true, default = 1.0)),
	α4 = $(@bind α4 Slider(0.0:0.01:10.0, show_value=true, default = 10.0)),\
	A = $(@bind A Slider(0.0:50.0, show_value=true, default = 0.0)),
	"""
end

# ╔═╡ cb3bf3b5-2f4f-4307-9bf3-82022d915d3c
Notch_model = @reaction_network begin
  (A,1.0),              R ↔ NR                   # NICD binds RBPJ
  (k1, k2),            M + R ↔ MR 		        # MITF binds RBPJ
# 	The epigenetic core regulation for the histone state of mir-222
	d, 		H4  + KDM5A  --> H0  + KDM5A        # Demethylation of active mark
	m,		H0  + PRC2   --> H27 + PRC2         # Methylation to get repressive mark
	1, 		H27 + KDM6A  --> H0  + KDM6A        # Demethylation of repressive mark
	1, 		H0  + KMT    --> H4  + KMT          # Methylation to get active mark
# 	Epigenetic feeback
	k0,			MR  --> MR  + KDM5A    				# MITF-RBPJ recruit KDM5A
	p, 		    H27 --> H27 + PRC2    				# PRC2 is enhenced by H27 mark
	kk,	 		H4  --> H4  + KDM6A   				# KDM6A is enhenced by H4 mark
	pp, 		H4  --> H4  + KMT     				# KMT is enhenced by H4 mark
	k, 		    H27 --> H27 + KDM5A   				# KDM5A is enhenced by H27 mark
	δ,			(PRC2, KDM5A, KDM6A, KMT) --> ∅     # Degradation of histone reader and writers
# 	Basal production rates for epigenetic writers
	α1,          ∅ --> (KDM5A,PRC2,KMT,KDM6A)
	# α2,			 ∅ --> (KMT) 
	# α3,			 ∅ --> (PRC2) 
	# α4,			 ∅ --> (KDM6A) 
end k0 k1 k2 d m p k pp kk δ α1 A


# ╔═╡ 626927b6-af8f-43f5-b53d-db5ad2dab579
latexify(convert(ODESystem, Notch_model))

# ╔═╡ 7cf91626-eb96-4d42-8531-19a12f7d45e6
params(Notch_model)

# ╔═╡ 4655614b-7d62-44a3-8b21-8f47ef95b0b6
Notch_model

# ╔═╡ e6d194fe-ffb2-4494-9a64-6c93a552ce2d
calculate_jacobian(convert(ODESystem, Notch_model))

# ╔═╡ 93512d98-b900-466a-a628-38b2fb162c9a
begin
	ts, cb  = make_cb([ti,tf], length(params(Notch_model)), activation)
	# u0 = initialization_con(50.0, 50.0, 50.0)
	# u0 = [17.0, 6.0, 473.0,  27.0,  21.0,  13.7,  23.8,  37.0,  12.4,  12.0, 37.0] # repressive mark
	u0 = [6.6,    0.0,  6.6,  43.4,   0.2,   948.3,  0.7,  1997.3,  399.4, 1.7, 1.2]
	tspan = (0.0,tmax)
	model_p = [k0, k1, k2, d, m, p, k, pp, kk, δ, α1, A]
	prob = ODEProblem(Notch_model, u0, tspan, model_p)
	sol = DifferentialEquations.solve(prob, callback=cb, tstops=ts)
	plot(sol, vars = [:MR, :KDM5A, :H4, :H27], lw  = 1.5)
end

# ╔═╡ 1fed2741-c239-45e1-bcb7-6f5722a92734
	@show u0

# ╔═╡ ee53e7c4-23f0-4cb0-85b6-6c39173a0c5e
md"""
## Test Bifurcation
"""

# ╔═╡ c4aa0019-6f78-4208-bdc3-b3cb7225682b
# begin
# 	odefun = ODEFunction(convert(ODESystem,Notch_model),jac=true)
# 	F = (u,p) -> odefun(u,p,0)      
# 	J = (u,p) -> odefun.jac(u,p,0)
# end

# ╔═╡ d0865215-e0d2-4149-bde9-59535bb60dae
begin
	p_idx = 4            # The index of the bifurcation parameter.
	p_span = (0.001,0.01)   # The parameter range for the bifurcation diagram.
	plot_var_idx = 5     # The index of the variable we plot in the bifurcation diagram.
end

# ╔═╡ 8aaa02aa-323a-45f2-bb35-84314ff141dd
# using BifurcationKit, LinearAlgebra, Setfield

# ╔═╡ 9956e025-49a6-471a-9a1e-9ab96236086b
# begin
# 	opts = ContinuationPar( dsmax = 0.05,        # Maximum arclength value of the pseudo-arc length continuation method.
# 	                        dsmin = 1e-4,        # Minimum arclength value of the pseudo-arc length continuation method.
# 	                        ds=0.001,            # Initial arclength value of the pseudo-arc length continuation method (should be positive).
# 	                        maxSteps = 100000,   # The maximum number of steps.
# 	                        pMin = p_span[1],    # Minimum p-vale (if hit, the method stops).
# 	                        pMax = p_span[2],    # Maximum p-vale (if hit, the method stops).
# 	                        detectBifurcation=3, # Value in {0,1,2,3} determening to what extent bofurcation points are detected (0 means nothing is done, 3 both them and there localisation are detected).
# 	                        newtonOptions = NewtonPar(tol = 1e-9, verbose = false, maxIter = 15)) #Parameters to the newton solver (when finding fixed points) see BifurcationKit documentation.
	                        
# 	DO = DeflationOperator( 2.0,    # Algorithm parameter required when using deflated continuation, see BifurcationKit documentation.
# 	                        dot,    # Algorithm parameter required when using deflated continuation, see BifurcationKit documentation.
# 	                        1.,     # Algorithm parameter required when using deflated continuation, see BifurcationKit documentation.
# 	                        [fill(0.,length(Notch_model.states))]); # Guess(es) of the fixed point for the initial parameter set. Do not need to be exact.
# end

# ╔═╡ 86b4effc-d714-4c03-93b6-eff4fed564d8
# params_input = setindex!(copy(model_p),p_span[1],p_idx)     # The input parameter values have to start at the first index of our parameter span.

# ╔═╡ e1a80cbe-4cd0-47c3-940e-488ffbe4af35
# begin
# 	branches, = continuation(F, J, params_input, (@lens _[p_idx]) ,opts , DO,             # Gives our input.
# 	    verbosity = 0, showplot=false,                                                    # We do not want to display, or plot, intermediary results.
# 	    printSolution=(x, p) -> x[plot_var_idx],                                          # How we wish to print the output in the diagram. Here we simply want the value of the target varriable.
# 	    perturbSolution = (x,p,id) -> (x  .+ 0.8 .* rand(length(x))),                     # Parameter for the continuation method, see BifurcationKit documentation.
# 	    callbackN = (x, f, J, res, iteration, itlinear, options; kwargs...) -> res <1e7)  # Parameter for the continuation method, see BifurcationKit documentation.
	
# end

# ╔═╡ 308dec25-3a51-49ac-9f39-0a248631f57c
md"""
Now let's create a function that takes in a parameter set and calculates the maximum of the predator population and the average of the prey population for those parameter values. To do this, we will make use of the remake function which creates a new ODEProblem, and use the p keyword argument to set the new parameters:
"""

# ╔═╡ 434c798d-b4af-47e2-b2e9-99d995a3a651
t = collect(range(0, stop=100, length=20000))

# ╔═╡ 42a6e0eb-9f73-4e7f-ba1e-02cccfb6416b
f1 = function (p)
  prob = remake(prob;p=p)
  sol = solve(prob,Tsit5();saveat=t)
  [mean(sol[6,:]), maximum(sol[6,:])]
end

# ╔═╡ c54de595-0526-41f6-97c3-e9fc6e247484


# ╔═╡ Cell order:
# ╟─347867d4-726a-4e83-977a-fd5e77512422
# ╠═d72bc2d4-b13f-450e-a24d-a056d7b48199
# ╠═dfb1f4ba-7ded-4b87-ac7d-4974349ed12d
# ╠═42fb32c8-10a6-4e75-940b-07193d4e549b
# ╠═8ea69aeb-13ef-4570-9150-4e23219fe8d4
# ╠═b5302649-d9be-47e3-9a11-d39370bb678a
# ╠═968fa436-802a-4ff8-bc9a-f225ac439e9d
# ╠═9e4cd586-71c0-4e7e-b958-1efe64f228b8
# ╠═2cc71eb3-df4f-4cf9-9682-1938959d58a2
# ╟─d8346212-0526-4423-be21-74cde773322d
# ╟─528975db-a04c-4e2c-8aee-be2083ea6056
# ╠═825f90df-6828-440e-b7b5-07106c1ae788
# ╠═36a50cf4-2b06-4e4b-9e66-469753e4f3b4
# ╟─8732d482-987a-40a8-9da7-8ec9523931bb
# ╠═cb3bf3b5-2f4f-4307-9bf3-82022d915d3c
# ╟─e7d396f4-ecf6-4d33-80ff-2fac41990906
# ╠═626927b6-af8f-43f5-b53d-db5ad2dab579
# ╟─7cf91626-eb96-4d42-8531-19a12f7d45e6
# ╟─35028063-80b0-4ab8-96b4-bc84af7389b9
# ╟─e86a64f0-1997-4149-b87a-e55b909cc153
# ╠═5bb6836f-698f-4276-a4bd-84c3260febfb
# ╟─a44a95ac-7cb7-4b4b-b64d-7dd43ec09629
# ╠═b5078a8f-8a6c-4257-ba4a-0f3ca0fac4c0
# ╠═4655614b-7d62-44a3-8b21-8f47ef95b0b6
# ╠═e6d194fe-ffb2-4494-9a64-6c93a552ce2d
# ╟─1eeae429-3c92-4bbf-9856-c7db46ccb95b
# ╠═3fef48e3-a6f9-4c64-bdea-41a750c4e6bf
# ╠═93512d98-b900-466a-a628-38b2fb162c9a
# ╠═1fed2741-c239-45e1-bcb7-6f5722a92734
# ╠═d47a550f-473e-47a7-bd89-e273f8f0854f
# ╠═a360a666-653d-4aed-b08d-88bdb955773a
# ╠═2b189db9-708a-49a1-8ef1-cc4a360fa8db
# ╟─ee53e7c4-23f0-4cb0-85b6-6c39173a0c5e
# ╠═c4aa0019-6f78-4208-bdc3-b3cb7225682b
# ╠═d0865215-e0d2-4149-bde9-59535bb60dae
# ╟─8aaa02aa-323a-45f2-bb35-84314ff141dd
# ╟─9956e025-49a6-471a-9a1e-9ab96236086b
# ╟─86b4effc-d714-4c03-93b6-eff4fed564d8
# ╟─e1a80cbe-4cd0-47c3-940e-488ffbe4af35
# ╟─308dec25-3a51-49ac-9f39-0a248631f57c
# ╠═42a6e0eb-9f73-4e7f-ba1e-02cccfb6416b
# ╠═3d428d36-cdf9-4f2b-99ee-14e7cbc1f39b
# ╠═434c798d-b4af-47e2-b2e9-99d995a3a651
# ╠═c54de595-0526-41f6-97c3-e9fc6e247484
