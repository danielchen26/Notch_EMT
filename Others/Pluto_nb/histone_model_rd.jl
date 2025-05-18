### A Pluto.jl notebook ###
# v0.16.0

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


# ╔═╡ e7d396f4-ecf6-4d33-80ff-2fac41990906
md"""
### From CRN to ODE
"""

# ╔═╡ 626927b6-af8f-43f5-b53d-db5ad2dab579
latexify(convert(ODESystem, Notch_model))

# ╔═╡ 7cf91626-eb96-4d42-8531-19a12f7d45e6
params(Notch_model)

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

# ╔═╡ 4655614b-7d62-44a3-8b21-8f47ef95b0b6
Notch_model

# ╔═╡ e6d194fe-ffb2-4494-9a64-6c93a552ce2d
calculate_jacobian(convert(ODESystem, Notch_model))

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


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Catalyst = "479239e8-5488-4da2-87a7-35f2df7eef83"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DiffEqBiological = "eb300fae-53e8-50a0-950c-e21f52c2b7e0"
DiffEqSensitivity = "41bf760c-e81c-5289-8e54-58b1f1f8abe2"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Catalyst = "~5.1.1"
DataFrames = "~0.21.8"
DiffEqBiological = "~4.3.0"
DiffEqSensitivity = "~6.35.0"
DifferentialEquations = "~6.18.0"
Images = "~0.24.1"
Latexify = "~0.13.5"
Plots = "~1.6.8"
PlutoUI = "~0.7.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractAlgebra]]
deps = ["InteractiveUtils", "LinearAlgebra", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "8a53c65dbc0423738c855abc42ea0c925ea1608b"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.11.2"

[[AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "345a14764e43fe927d6f5c250fe4c8e4664e6ee8"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "2.4.0"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[ArrayInterface]]
deps = ["LinearAlgebra", "Requires", "SparseArrays"]
git-tree-sha1 = "a2a1884863704e0414f6f164a1f6f4a2a62faf4e"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "2.14.17"

[[ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "951c3fc1ff93497c88fb1dfa893f4de55d0b38e3"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.3.8"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "a4d07a1c313392a77042855df46c5f534076fab9"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.0"

[[AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "d127d5e4d86c7680b20c35d40b503c74b9a39b5e"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.4"

[[BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "eaf98fa821ab26c5825fa3a054db735755c335da"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.15.15"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BinaryProvider]]
deps = ["Libdl", "Logging", "SHA"]
git-tree-sha1 = "ecdec412a9abc8db54c0efc5548c64dfce072058"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.10"

[[BoundaryValueDiffEq]]
deps = ["BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NLsolve", "Reexport", "SparseArrays"]
git-tree-sha1 = "fe34902ac0c3a35d016617ab7032742865756d7d"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "2.7.1"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f2202b55d816427cd385a9a4f3ffb226bee80f99"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+0"

[[CanonicalTraits]]
deps = ["MLStyle"]
git-tree-sha1 = "f959d0e7164fb0262b02abecb93cf42b9a9f3188"
uuid = "a603d957-0e48-4f86-8fbd-0b7bc66df689"
version = "0.2.4"

[[CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[Catalyst]]
deps = ["Catlab", "Latexify", "MacroTools", "ModelingToolkit", "Reexport"]
git-tree-sha1 = "b50f703b5835f6d277d0096ba0e3bc0c52caf8ae"
uuid = "479239e8-5488-4da2-87a7-35f2df7eef83"
version = "5.1.1"

[[CategoricalArrays]]
deps = ["DataAPI", "Future", "JSON", "Missings", "Printf", "Statistics", "StructTypes", "Unicode"]
git-tree-sha1 = "2ac27f59196a68070e132b25713f9a5bbc5fa0d2"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.8.3"

[[Catlab]]
deps = ["AutoHashEquals", "Compat", "Compose", "DataStructures", "FunctionWrappers", "GeneralizedGenerated", "JSON", "LightGraphs", "LightXML", "LinearAlgebra", "Logging", "MLStyle", "Parameters", "Pkg", "PrettyTables", "Random", "Reexport", "Requires", "SparseArrays", "StaticArrays", "Statistics", "StructArrays"]
git-tree-sha1 = "4503d0ec30bdea6ea1a49ad6e4cbf8e6a38be5cf"
uuid = "134e5e36-593f-5add-ad60-77f754baafbe"
version = "0.9.2"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4b28f88cecf5d9a07c85b9ce5209a361ecaff34a"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "0.9.45"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "9995eb3977fbf67b86d0a0a0508e83017ded03f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.14.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "32a2b8af383f11cbb65803883837a149d10dfe8a"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.10.12"

[[ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "a66a8e024807c4b3d186eb1cab2aff3505271f8e"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.6"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "4866e381721b30fac8dda4c8cb1d9db45c8d2994"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.37.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "c6461fc7c35a4bb8d00905df7adafcff1fe3a6bc"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.2"

[[ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[ConsoleProgressMonitor]]
deps = ["Logging", "ProgressMeter"]
git-tree-sha1 = "3ab7b2136722890b9af903859afcf457fa3059e8"
uuid = "88cd18e8-d9cc-4ea6-8889-5259c0d15c8b"
version = "0.1.2"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "6d1c23e740a586955645500bbec662476204a52c"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.1"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataFrames]]
deps = ["CategoricalArrays", "Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "Missings", "PooledArrays", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ecd850f3d2b815431104252575e7307256121548"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "0.21.8"

[[DataStructures]]
deps = ["InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "88d48e133e6d3dd68183309877eac74393daa7eb"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.17.20"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelayDiffEq]]
deps = ["DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "Roots", "UnPack"]
git-tree-sha1 = "ede535e084589ca2ed626f4179337ac7fd610e7f"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.24.1"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "ConsoleProgressMonitor", "DataStructures", "Distributed", "DocStringExtensions", "FunctionWrappers", "IterativeSolvers", "IteratorInterfaceExtensions", "LabelledArrays", "LinearAlgebra", "Logging", "LoggingExtras", "MuladdMacro", "Parameters", "Printf", "ProgressLogging", "RecipesBase", "RecursiveArrayTools", "RecursiveFactorization", "Requires", "Roots", "SparseArrays", "StaticArrays", "Statistics", "SuiteSparse", "TableTraits", "TerminalLoggers", "TreeViews", "ZygoteRules"]
git-tree-sha1 = "4dde9a142bc3780f0216e564673f477ec6f6a6df"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.44.3"

[[DiffEqBiological]]
deps = ["Compat", "DataStructures", "DiffEqBase", "DiffEqJump", "DynamicPolynomials", "HomotopyContinuation", "Latexify", "LinearAlgebra", "MacroTools", "Parameters", "RecipesBase", "Reexport", "SparseArrays", "SymEngine"]
git-tree-sha1 = "b28f4547b13ab348f0de7ec9eb369cbf6f72fa60"
uuid = "eb300fae-53e8-50a0-950c-e21f52c2b7e0"
version = "4.3.0"

[[DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "NLsolve", "OrdinaryDiffEq", "RecipesBase", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "63afe1bc146559c7329b6429713916c7f9e61d2f"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.13.5"

[[DiffEqFinancial]]
deps = ["DiffEqBase", "DiffEqNoiseProcess", "LinearAlgebra", "Markdown", "RandomNumbers"]
git-tree-sha1 = "db08e0def560f204167c58fd0637298e13f58f73"
uuid = "5a0ffddc-d203-54b0-88ba-2c03c0fc2e67"
version = "2.4.0"

[[DiffEqJump]]
deps = ["ArrayInterface", "Compat", "DataStructures", "DiffEqBase", "FunctionWrappers", "LinearAlgebra", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "97c44c5f94827bb4605e4c72b652629a9d6d10d4"
uuid = "c894b116-72e5-5b58-be3c-e6d8d4ac2b12"
version = "6.10.0"

[[DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "LinearAlgebra", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArrays", "Statistics"]
git-tree-sha1 = "d6839a44a268c69ef0ed927b22a6f43c8a4c2e73"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.9.0"

[[DiffEqPhysics]]
deps = ["DiffEqBase", "DiffEqCallbacks", "ForwardDiff", "LinearAlgebra", "Printf", "Random", "RecipesBase", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "8f23c6f36f6a6eb2cbd6950e28ec7c4b99d0e4c9"
uuid = "055956cb-9e8b-5191-98cc-73ae4a59e68a"
version = "3.9.0"

[[DiffEqSensitivity]]
deps = ["Adapt", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "Distributions", "FFTW", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Parameters", "QuadGK", "QuasiMonteCarlo", "RecursiveArrayTools", "Requires", "ReverseDiff", "SharedArrays", "Statistics", "StochasticDiffEq", "Tracker", "Zygote", "ZygoteRules"]
git-tree-sha1 = "ba439d161a82a50dbf9f3297281ecbfa9a4c1642"
uuid = "41bf760c-e81c-5289-8e54-58b1f1f8abe2"
version = "6.35.0"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "7220bc21c33e990c14f4a9a319b1d242ebc5b269"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.3.1"

[[DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqFinancial", "DiffEqJump", "DiffEqNoiseProcess", "DiffEqPhysics", "DimensionalPlotRecipes", "LinearAlgebra", "MultiScaleArrays", "OrdinaryDiffEq", "ParameterizedFunctions", "Random", "RecursiveArrayTools", "Reexport", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "ececc535bd2aa55a520131d955639288704e3851"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "6.18.0"

[[DimensionalPlotRecipes]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "af883a26bbe6e3f5f778cb4e1b81578b534c32a6"
uuid = "c619ae07-58cd-5f6d-b883-8f17bd6a98f9"
version = "1.2.0"

[[Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "9f46deb4d4ee4494ffb5a40a27a2aced67bdd838"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.4"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "9c41285c57c6e0d73a21ed4b65f6eec34805f937"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.23.8"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[DoubleFloats]]
deps = ["GenericLinearAlgebra", "LinearAlgebra", "Polynomials", "Printf", "Quadmath", "Random", "Requires", "SpecialFunctions"]
git-tree-sha1 = "1c962cf7e75c09a5f1fbf504df7d6a06447a1129"
uuid = "497a8b3b-efae-58df-a0af-a86822472b78"
version = "1.1.23"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "05b68e727a192783be0b34bd8fee8f678505c0bf"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.3.20"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[ElasticArrays]]
deps = ["Adapt"]
git-tree-sha1 = "a0fcc1bb3c9ceaf07e1d0529c9806ce94be6adf9"
uuid = "fdbdab4c-e67f-52f5-8c3f-e7b388dad3d4"
version = "1.2.9"

[[EllipsisNotation]]
git-tree-sha1 = "18ee049accec8763be17a933737c1dd0fdf8673a"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.0.0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[ExponentialUtilities]]
deps = ["LinearAlgebra", "Printf", "Requires", "SparseArrays"]
git-tree-sha1 = "3546dbb394f7b95d805d6b1d9ccf63310ba59557"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.8.2"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "70a0cfd9b1c86b0209e38fbfe6d8231fd606eeaf"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.1"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "463cb335fa22c4ebacfd1faba5fde14edb80d96c"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.5"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "3c041d2ac0a52a12a27af2782b34900d9c3ee68c"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.11.1"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "4863cbb7910079369e258dee4add9d06ead5063a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.8.14"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[FixedPolynomials]]
deps = ["LinearAlgebra", "MultivariatePolynomials", "Test"]
git-tree-sha1 = "32ec33fe4387eb165461b7141d4be732097086ad"
uuid = "3dd14ad9-0029-526e-86e9-8aa0bdd2ab0d"
version = "0.4.0"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "b5e930ac60b613ef3406da6d4f42c35d8dc51419"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.19"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GMP_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "781609d7-10c4-51f6-84f2-b8444358ff6d"

[[GR]]
deps = ["Base64", "DelimitedFiles", "HTTP", "JSON", "LinearAlgebra", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "cd0f34bd097d4d5eb6bbe01778cf8a7ed35f29d9"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.52.0"

[[GeneralizedGenerated]]
deps = ["CanonicalTraits", "DataStructures", "JuliaVariables", "MLStyle"]
git-tree-sha1 = "7dd404baf79b28f117917633f0cc1d2976c1fd9f"
uuid = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
version = "0.2.8"

[[GenericLinearAlgebra]]
deps = ["LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "eddbb6ee8fe2c3244a2c973874a3179c3c4d3ac5"
uuid = "14197337-ba66-59df-a3e3-ca00e7dcff7a"
version = "0.2.6"

[[GenericSVD]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "62909c3eda8a25b5673a367d1ad2392ebb265211"
uuid = "01680d73-4ee2-5a08-a1aa-533608c188bb"
version = "0.3.0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "82853ebc70db4f5a3084853738c68fd497b22c7c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.3.10"

[[GeometryTypes]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "07194161fe4e181c6bf51ef2e329ec4e7d050fc4"
uuid = "4d00f742-c7ba-57c2-abde-4428a4b178cb"
version = "0.8.4"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7bf67e9a481712b3dbe9cb3dac852dc4b1162e02"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+0"

[[Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "2c1cf4df419938ece72de17f368a021ee162762e"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.0"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "MbedTLS", "Sockets"]
git-tree-sha1 = "c7ec02c4c6a039a98a15f955462cd7aea5df4508"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.8.19"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "8a954fed8ac097d5be04921d595f741115c1b2ad"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+0"

[[HomotopyContinuation]]
deps = ["DoubleFloats", "DynamicPolynomials", "ElasticArrays", "FixedPolynomials", "Latexify", "LinearAlgebra", "MixedSubdivisions", "MultivariatePolynomials", "OrderedCollections", "Parameters", "PrettyTables", "Printf", "ProjectiveVectors", "Random", "StaticArrays", "StaticPolynomials", "Test", "TreeViews"]
git-tree-sha1 = "794f8b33b2bbaa5c085a4bef5f684cc026b3ea91"
uuid = "f213a82b-91d6-5c5d-acf7-10f1c761b327"
version = "1.4.4"

[[Hwloc]]
deps = ["Hwloc_jll"]
git-tree-sha1 = "ffdcd4272a7cc36442007bca41aa07ca3cc5fda4"
uuid = "0e44f5e4-bd66-52a0-8798-143a42290a1d"
version = "1.3.0"

[[Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3395d4d4aeb3c9d31f5929d32760d8baeee88aaf"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.5.0+0"

[[IRTools]]
deps = ["InteractiveUtils", "MacroTools", "Test"]
git-tree-sha1 = "90ee39f9beaaa186e4968417ea2b8ed5673c91c0"
uuid = "7869d1d1-7146-5819-86e3-90919afe41df"
version = "0.3.3"

[[IdentityRanges]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be8fcd695c4da16a1d6d0cd213cb88090a150e3b"
uuid = "bbac6d45-d8f3-5730-bfe4-7a449cd117ca"
version = "0.3.1"

[[IfElse]]
git-tree-sha1 = "28e837ff3e7a6c3cdb252ce49fb412c8eb3caeef"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.0"

[[ImageAxes]]
deps = ["AxisArrays", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "794ad1d922c432082bc1aaa9fa8ffbd1fe74e621"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.9"

[[ImageContrastAdjustment]]
deps = ["ColorVectorSpace", "ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "2e6084db6cccab11fe0bc3e4130bd3d117092ed9"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.7"

[[ImageCore]]
deps = ["AbstractFFTs", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "db645f20b59f060d8cfae696bc9538d13fd86416"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.8.22"

[[ImageDistances]]
deps = ["ColorVectorSpace", "Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "6378c34a3c3a216235210d19b9f495ecfff2f85f"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.13"

[[ImageFiltering]]
deps = ["CatIndices", "ColorVectorSpace", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageCore", "LinearAlgebra", "OffsetArrays", "Requires", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "bf96839133212d3eff4a1c3a80c57abc7cfbf0ce"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.6.21"

[[ImageIO]]
deps = ["FileIO", "Netpbm", "OpenEXR", "PNGFiles", "TiffImages", "UUIDs"]
git-tree-sha1 = "13c826abd23931d909e4c5538643d9691f62a617"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.5.8"

[[ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[ImageMagick_jll]]
deps = ["JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1c0a2295cca535fabaf2029062912591e9b61987"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.10-12+3"

[[ImageMetadata]]
deps = ["AxisArrays", "ColorVectorSpace", "ImageAxes", "ImageCore", "IndirectArrays"]
git-tree-sha1 = "ae76038347dc4edcdb06b541595268fca65b6a42"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.5"

[[ImageMorphology]]
deps = ["ColorVectorSpace", "ImageCore", "LinearAlgebra", "TiledIteration"]
git-tree-sha1 = "68e7cbcd7dfaa3c2f74b0a8ab3066f5de8f2b71d"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.2.11"

[[ImageQualityIndexes]]
deps = ["ColorVectorSpace", "ImageCore", "ImageDistances", "ImageFiltering", "OffsetArrays", "Statistics"]
git-tree-sha1 = "1198f85fa2481a3bb94bf937495ba1916f12b533"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.2.2"

[[ImageShow]]
deps = ["Base64", "FileIO", "ImageCore", "OffsetArrays", "Requires", "StackViews"]
git-tree-sha1 = "832abfd709fa436a562db47fd8e81377f72b01f9"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.1"

[[ImageTransformations]]
deps = ["AxisAlgorithms", "ColorVectorSpace", "CoordinateTransformations", "IdentityRanges", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "e4cc551e4295a5c96545bb3083058c24b78d4cf0"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.8.13"

[[Images]]
deps = ["AxisArrays", "Base64", "ColorVectorSpace", "FileIO", "Graphics", "ImageAxes", "ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageShow", "ImageTransformations", "IndirectArrays", "OffsetArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "8b714d5e11c91a0d945717430ec20f9251af4bd2"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.24.1"

[[Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[IndirectArrays]]
git-tree-sha1 = "c2a145a145dc03a7620af1444e0264ef907bd44f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "0.5.1"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["AxisAlgorithms", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "1e0e51692a3a77f1eeb51bf741bdd0439ed210e7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.2"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "704eee044a41b0e7e8417f7dd2a6b6b5361afd5f"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.8.5"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[JuliaVariables]]
deps = ["MLStyle", "NameResolution"]
git-tree-sha1 = "49fb3cb53362ddadb4415e9b73926d6b40709e70"
uuid = "b14d175d-62b4-44ba-8fb7-3064adc8c3ec"
version = "0.2.4"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[LabelledArrays]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "bdde43e002847c34c206735b1cf860bc3abd35e7"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.6.4"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "864527aa4d14c893fb8c51d48ef314410c88c7b9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.13.5"

[[LatinHypercubeSampling]]
deps = ["Random", "StableRNGs", "StatsBase", "Test"]
git-tree-sha1 = "42938ab65e9ed3c3029a8d2c58382ca75bdab243"
uuid = "a5e1c1ea-c99a-51d3-a14d-a9a37257b02d"
version = "1.8.0"

[[LatticeRules]]
deps = ["Random"]
git-tree-sha1 = "7f5b02258a3ca0221a6a9710b0a0a2e8fb4957fe"
uuid = "73f95e8e-ec14-4e6a-8b18-0d2e271c4e55"
version = "0.0.1"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "71be1eb5ad19cb4f61fa8c73395c0338fd092ae0"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.1.2"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

[[LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "e129d9391168c677cd4800f5c0abb1ed8cb3794f"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.0"

[[LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "7bd5f6565d80b6bf753738d2bc40a5dfea072070"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.2.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "dfeda1c1130990428720de0024d4516b1902ce98"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.7"

[[LoopVectorization]]
deps = ["ArrayInterface", "DocStringExtensions", "IfElse", "LinearAlgebra", "OffsetArrays", "Requires", "SLEEFPirates", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "5f275de503982d59bd82eb1e4fbc273f55a72dee"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.10.0"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[MLStyle]]
git-tree-sha1 = "594e189325f66e23a8818e5beb11c43bb0141bcd"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.10"

[[MPC_jll]]
deps = ["GMP_jll", "Libdl", "MPFR_jll", "Pkg"]
git-tree-sha1 = "583d9bc863ad491571369212b9a8047219a0015d"
uuid = "2ce0c516-f11f-5db3-98ad-e0e1048fbd70"
version = "1.1.0+0"

[[MPFR_jll]]
deps = ["Artifacts", "GMP_jll", "Libdl"]
uuid = "3a97d323-0669-5f0c-9066-3539efd106a3"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "5a5bc6bf062f0f95e62d0fe0a2d99699fed82dd9"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.8"

[[MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f8c673ccc215eb50fcadb285f522420e29e69e1c"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "0.4.5"

[[MixedSubdivisions]]
deps = ["LinearAlgebra", "MultivariatePolynomials", "ProgressMeter", "StaticArrays"]
git-tree-sha1 = "8671af51cd89bbca7fe6856634dcc250fb3038c9"
uuid = "291d046c-3347-11e9-1e74-c3d251d406c6"
version = "1.1.1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["ExprTools"]
git-tree-sha1 = "748f6e1e4de814b101911e64cc12d83a6af66782"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.2"

[[ModelingToolkit]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqJump", "DiffRules", "Distributed", "DocStringExtensions", "GeneralizedGenerated", "IfElse", "LabelledArrays", "Latexify", "Libdl", "LightGraphs", "LinearAlgebra", "MacroTools", "NaNMath", "RecursiveArrayTools", "Requires", "SafeTestsets", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews", "UnPack", "Unitful"]
git-tree-sha1 = "a2777851fd21be36882ddafe5b5e2036cdc14636"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "3.20.1"

[[MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[MultiScaleArrays]]
deps = ["DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "SparseDiffTools", "Statistics", "StochasticDiffEq", "TreeViews"]
git-tree-sha1 = "258f3be6770fe77be8870727ba9803e236c685b8"
uuid = "f9640e96-87f6-5992-9c3b-0743c6a49ffa"
version = "1.8.1"

[[MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "45c9940cec79dedcdccc73cc6dd09ea8b8ab142c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.3.18"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3927848ccebcc165952dc0d9ac9aa274a87bfe01"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.20"

[[NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "144bab5b1443545bc4e791536c9f1eacb4eed06a"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.1"

[[NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[NNlib]]
deps = ["BinaryProvider", "Libdl", "LinearAlgebra", "Requires", "Statistics"]
git-tree-sha1 = "d9f196d911f55aeaff11b11f681b135980783824"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.6.6"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NameResolution]]
deps = ["PrettyPrint"]
git-tree-sha1 = "1a0fa0e9613f46c9b8c11eee38ebb4f590013c5e"
uuid = "71a1bf82-56d0-4bbc-8a3c-48b961074391"
version = "0.1.5"

[[Netpbm]]
deps = ["ColorVectorSpace", "FileIO", "ImageCore"]
git-tree-sha1 = "09589171688f0039f13ebe0fdcc7288f50228b52"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.1"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "c0e9e582987d36d5a61e650e6e543b9e44d9914b"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.7"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Optim]]
deps = ["Compat", "FillArrays", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "7863df65dbb2a0fa8f85fcaf0a41167640d2ebed"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.4.1"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "ExponentialUtilities", "FastClosures", "FiniteDiff", "ForwardDiff", "GenericSVD", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "NLsolve", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "e6806009a3a5e74c5a6d44365c15b967daa7f798"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "5.42.3"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "95a4038d1011dfdbde7cecd2ad0ac411e53ab1bc"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.10.1"

[[PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "e14c485f6beee0c7a8dcf6128bf70b85f1fe201e"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.9"

[[PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "646eed6f6a5d8df6708f15ea7e02a7a2c4fe4800"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.10"

[[ParameterizedFunctions]]
deps = ["DataStructures", "DiffEqBase", "Latexify", "LinearAlgebra", "ModelingToolkit", "Reexport"]
git-tree-sha1 = "291279c720121d7f5c6a145726bb94da79e9b42c"
uuid = "65888b18-ceab-5e60-b2b9-181511a3b968"
version = "5.6.0"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "2276ac65f1e236e0a6ea70baff3f62ad4c625345"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.2"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "2537ed3c0ed5e03896927187f5f2ee6a4ab342db"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.14"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "GeometryTypes", "JSON", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "c499e18bbeab024f6de0e0ae285554d153eeb5c5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.6.8"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "Logging", "Markdown", "Random", "Suppressor"]
git-tree-sha1 = "45ce174d36d3931cd4e37a47f93e07d1455f038d"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.1"

[[PoissonRandom]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "44d018211a56626288b5d3f8c6497d28c26dc850"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.0"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "0bbfdcd8cda81b8144de4be8a67f5717e959a005"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "2.0.14"

[[PooledArrays]]
deps = ["DataAPI"]
git-tree-sha1 = "b1333d4eced1826e15adbdf01a4ecaccca9d353c"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "0.5.3"

[[PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[PrettyPrint]]
git-tree-sha1 = "632eb4abab3449ab30c5e1afaa874f0b98b586e4"
uuid = "8162dcfd-2161-5ef2-ae6c-7681170c5f98"
version = "0.2.0"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Parameters", "Reexport", "Tables"]
git-tree-sha1 = "8458dc04a493ae5c2fed3796c1d3117972c69694"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "0.9.1"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "afadeba63d90ff223a6a48d2009434ecee2ec9e8"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.1"

[[ProjectiveVectors]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4800d531755d189b1a28c0eb48bc9d9311273a70"
uuid = "01f381cc-face-5a0a-ade9-ef63dc65d628"
version = "1.1.3"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[Quadmath]]
deps = ["Printf", "Random", "Requires"]
git-tree-sha1 = "5a8f74af8eae654086a1d058b4ec94ff192e3de0"
uuid = "be4d8f0f-7fa4-5f49-b795-2f01399ab2dd"
version = "0.5.5"

[[QuasiMonteCarlo]]
deps = ["Distributions", "LatinHypercubeSampling", "LatticeRules", "Sobol"]
git-tree-sha1 = "e1ffcdd1f7a81db4376b474a441691084d5b7a82"
uuid = "8a4e6c94-4038-4cdc-81c3-7e6ffdb2a71b"
version = "0.2.3"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Random123]]
deps = ["Libdl", "Random", "RandomNumbers"]
git-tree-sha1 = "0e8b146557ad1c6deb1367655e052276690e71a3"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.4.2"

[[RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[Ratios]]
deps = ["Requires"]
git-tree-sha1 = "01d341f502250e81f6fec0afe662aa861392a3aa"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.2"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "4a325c9bcc2d8e62a8f975b9666d0251d53b63b9"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.1.13"

[[RecursiveArrayTools]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "b3f4e34548b3d3d00e5571fd7bc0a33980f01571"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.11.4"

[[RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization"]
git-tree-sha1 = "2e1a88c083ebe8ba69bc0b0084d4b4ba4aa35ae0"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.1.13"

[[Reexport]]
deps = ["Pkg"]
git-tree-sha1 = "7b1d07f411bc8ddb7977ec7f377b97b158514fe0"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "0.2.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[ReverseDiff]]
deps = ["DiffResults", "DiffRules", "ForwardDiff", "FunctionWrappers", "LinearAlgebra", "MacroTools", "NaNMath", "Random", "SpecialFunctions", "StaticArrays", "Statistics"]
git-tree-sha1 = "63ee24ea0689157a1113dbdab10c6cb011d519c4"
uuid = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
version = "1.9.0"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[Roots]]
deps = ["CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "e8102d671f2343d2e9f12af0ecbd85fd5f86c4fa"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "1.3.4"

[[Rotations]]
deps = ["LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "2ed8d8a16d703f900168822d83699b8c3c1a5cd8"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.0.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SLEEFPirates]]
deps = ["IfElse", "Libdl", "VectorizationBase"]
git-tree-sha1 = "3d44bb7517298fd262915924fdc1645c61a6ef17"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.8"

[[SafeTestsets]]
deps = ["Test"]
git-tree-sha1 = "36ebc5622c82eb9324005cc75e7e2cc51181d181"
uuid = "1bc83da4-3b8d-516f-aca4-4fe02f6d838f"
version = "0.0.1"

[[SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "0b5b04995a3208f82ae69a17389c60b762000793"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.13.6"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "fca29e68c5062722b5b4435594c3d1ba557072a3"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.7.1"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "ee010d8f103468309b8afac4abb9be2e18ff1182"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "0.3.2"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sobol]]
deps = ["DelimitedFiles", "Random"]
git-tree-sha1 = "5a74ac22a9daef23705f010f72c81d6925b19df8"
uuid = "ed01d8cd-4d21-5b2a-85b4-cc3bdc58bad4"
version = "1.5.0"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures", "Random", "Test"]
git-tree-sha1 = "03f5898c9959f8115e30bc7226ada7d0df554ddd"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "0.3.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "LightGraphs", "LinearAlgebra", "Requires", "SparseArrays", "VertexSafeGraphs"]
git-tree-sha1 = "be20320958ccd298c98312137a5ebe75a654ebc8"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.13.2"

[[SpecialFunctions]]
deps = ["OpenSpecFun_jll"]
git-tree-sha1 = "d8d8b8a9f4119829410ecd706da4cc8594a1e020"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "0.10.3"

[[StableRNGs]]
deps = ["Random", "Test"]
git-tree-sha1 = "3be7d49667040add7ee151fefaf1f8c04c8c8276"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.0"

[[StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "da4cf579416c81994afd6322365d00916c79b8ae"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "0.12.5"

[[StaticPolynomials]]
deps = ["LinearAlgebra", "MultivariatePolynomials", "StaticArrays"]
git-tree-sha1 = "0a30b1862911d0889494c1b8b3d6f028860f5577"
uuid = "62e018b1-6e46-5407-a5a7-97d4fbcae734"
version = "1.3.4"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8cbbc098554648c84f79a463c9ff0fd277144b6c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.10"

[[StatsFuns]]
deps = ["LogExpFunctions", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "30cd8c360c54081f806b1ee14d2eecbef3c04c49"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.8"

[[SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "3df66a4a9ba477bea5cb10a3ec732bb48a2fc27d"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.6.4"

[[StochasticDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqJump", "DiffEqNoiseProcess", "FillArrays", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "ea25d14313775646f9aa4b2158156fb1b3815001"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.25.0"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "Tables"]
git-tree-sha1 = "8099ed9fb90b6e754d6ba8c6ed8670f010eadca0"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.4.4"

[[StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "8445bf99a36d703a09c601f9a57e2f83000ef2ae"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.7.3"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"

[[Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "Reexport", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "0fc68b00d45b3882174d21a5841a1049769b9fd8"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.2.6"

[[Sundials_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl", "OpenBLAS_jll", "Pkg", "SuiteSparse_jll"]
git-tree-sha1 = "013ff4504fc1d475aa80c63b455b6b3a58767db2"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.0+1"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[SymEngine]]
deps = ["Compat", "Libdl", "LinearAlgebra", "RecipesBase", "SpecialFunctions", "SymEngine_jll"]
git-tree-sha1 = "c2d81981a69ea69713747c0ce9e0f95883d3ff9b"
uuid = "123dc426-2d89-5057-bbad-38513e3affd8"
version = "0.8.5"

[[SymEngine_jll]]
deps = ["GMP_jll", "Libdl", "MPC_jll", "MPFR_jll", "Pkg"]
git-tree-sha1 = "4dacada8e05ac49eb768219f8d02bc6b608627fb"
uuid = "3428059b-622b-5399-b16f-d347a77089a4"
version = "0.6.0+1"

[[SymbolicUtils]]
deps = ["AbstractAlgebra", "Combinatorics", "DataStructures", "NaNMath", "SpecialFunctions", "TimerOutputs"]
git-tree-sha1 = "cd230ab5f02844155415aad28e8474fe459fe366"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.5.2"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "1162ce4a6c4b7e31e0e6b14486a6986951c73be9"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.2"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[TerminalLoggers]]
deps = ["LeftChildRightSiblingTrees", "Logging", "Markdown", "Printf", "ProgressLogging", "UUIDs"]
git-tree-sha1 = "d620a061cb2a56930b52bdf5cf908a5c4fa8e76a"
uuid = "5d786b92-1e48-4d6f-9151-6b4477ca9bed"
version = "0.1.4"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[ThreadingUtilities]]
deps = ["VectorizationBase"]
git-tree-sha1 = "64a2be7c73951d7c402eb40a16055e5e6fdda468"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.2.3"

[[TiffImages]]
deps = ["ColorTypes", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "OrderedCollections", "PkgVersion", "ProgressMeter"]
git-tree-sha1 = "632a8d4dbbad6627a4d2d21b1c6ebcaeebb1e1ed"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.4.2"

[[TiledIteration]]
deps = ["OffsetArrays"]
git-tree-sha1 = "52c5f816857bfb3291c7d25420b1f4aca0a74d18"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.3.0"

[[TimeZones]]
deps = ["Dates", "Future", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "6c9040665b2da00d30143261aea22c7427aada1c"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.5.7"

[[TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "209a8326c4f955e2442c07b56029e88bb48299c7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.12"

[[Tracker]]
deps = ["Adapt", "DiffRules", "ForwardDiff", "LinearAlgebra", "MacroTools", "NNlib", "NaNMath", "Printf", "Random", "Requires", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "bf4adf36062afc921f251af4db58f06235504eff"
uuid = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
version = "0.2.16"

[[TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a981a8ef8714cba2fd9780b22fd7a469e7aaf56d"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.9.0"

[[VectorizationBase]]
deps = ["ArrayInterface", "Hwloc", "IfElse", "Libdl", "LinearAlgebra"]
git-tree-sha1 = "9f27ddc74f04319574749bb0a5b3f231a5e7db16"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.16.2"

[[VertexSafeGraphs]]
deps = ["LightGraphs"]
git-tree-sha1 = "b9b450c99a3ca1cc1c6836f560d8d887bcbe356e"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.1.2"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "59e2ad8fd1591ea019a5259bd012d7aee15f995c"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.3"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[Zygote]]
deps = ["DiffRules", "FFTW", "FillArrays", "ForwardDiff", "IRTools", "InteractiveUtils", "LinearAlgebra", "MacroTools", "NNlib", "NaNMath", "Random", "Requires", "SpecialFunctions", "Statistics", "ZygoteRules"]
git-tree-sha1 = "f8329b595c465caf3ca87c4f744e6041a4983e43"
uuid = "e88e6eb3-aa80-5325-afca-941959d7151f"
version = "0.4.8"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "9e7a1e8ca60b742e508a315c17eef5211e7fbfd7"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.1"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
