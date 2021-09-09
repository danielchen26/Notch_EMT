### A Pluto.jl notebook ###
# v0.14.5

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

# ╔═╡ 42fb32c8-10a6-4e75-940b-07193d4e549b
using DiffEqBiological:steady_states

# ╔═╡ 347867d4-726a-4e83-977a-fd5e77512422
md"""
# Notch Epigenetics
"""

# ╔═╡ d72bc2d4-b13f-450e-a24d-a056d7b48199
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ 58c83dca-5e6c-4654-95ba-7e791fb93de6


# ╔═╡ b5302649-d9be-47e3-9a11-d39370bb678a
plotly()

# ╔═╡ d8346212-0526-4423-be21-74cde773322d
md"""
## Defining  the model
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

# ╔═╡ e456b730-ce1a-49d9-af9f-0818fd5b9d7c
begin
	# p  = params(Notch_model)
	md"""
	N = $(@bind N Slider(0.0:100.0, show_value=true, default = 1.0)) 
	k0 = $(@bind k0 Slider(0.0:100.0, show_value=true, default = 2.0))
	k1 = $(@bind k1 Slider(0.0:100.0, show_value=true, default = 20.0))
	k2 = $(@bind k2 Slider(0.0:100.0, show_value=true, default = 20.0))
	d1 = $(@bind d1 Slider(0.0:100.0, show_value=true, default = 20.0))
	d0 = $(@bind d0 Slider(0.0:100.0, show_value=true, default = 20.0))
	m0 = $(@bind m0 Slider(0.0:100.0, show_value=true, default = 20.0))
	m1 = $(@bind m1 Slider(0.0:100.0, show_value=true, default = 20.0))
	p = $(@bind p Slider(0.0:100.0, show_value=true, default = 20.0))
	k = $(@bind k Slider(0.0:100.0, show_value=true, default = 20.0))
	pp = $(@bind pp Slider(0.0:100.0, show_value=true, default = 20.0))
	kk = $(@bind kk Slider(0.0:100.0, show_value=true, default = 20.0))
	δ = $(@bind δ Slider(0.0:100.0, show_value=true, default = 1.0))
	v1 = $(@bind v1 Slider(0.0:100.0, show_value=true, default = 10.0))
	Kd1 = $(@bind Kd1 Slider(0.0:100.0, show_value=true, default = 50.0))
	v2 = $(@bind v2 Slider(0.0:100.0, show_value=true, default = 10.0))
	Kd2 = $(@bind Kd2 Slider(0.0:100.0, show_value=true, default = 50.0))
	
	"""
end

# ╔═╡ cb3bf3b5-2f4f-4307-9bf3-82022d915d3c
Notch_model = @reaction_network begin
	
# 	The competitve binding between NICD and MITF
# 	model NICD as an external input
	
  (N,k0),              R ↔ NR               # NICD binds RBPJ           
  (k1, k2),            M + R ↔ MR 		    # MITF binds RBPJ
   hill(M,v1,Kd1,2),   M --> M + R          # MITF form homodimer to activate RBPJ
   hillr(N,v2,Kd2,1),   N + M --> N         # NICD transcriptionally inhibit MITF
	
# 	The epigenetic core regulation for the histone state of mir-222
	
	d1, 		H4  + KDM5A  --> H0  + KDM5A  # Demethylation of active mark
	m0,			H0  + PRC2   --> H27 + PRC2   # Methylation to get repressive mark
	d0, 		H27 + KDM6A  --> H0  + KDM6A  # Demethylation of repressive mark
	m1, 		H0  + KMT    --> H4  + KMT    # Methylation to get active mark
	
# 	Epigenetic feeback 
	k0,			MR --> MR + KDM5A    				# MITF-RBPJ recruit KDM5A
	p, 		    H27 --> H27 + PRC2    				# PRC2 is enhenced by H27 mark     
	kk,	 		H4  --> H4  + KDM6A   				# KDM6A is enhenced by H4 mark
	pp, 		H4  --> H4  + KMT     				# KMT is enhenced by H4 mark
	k, 		    H27 --> H27 + KDM5A   				# KDM5A is enhenced by H27 mark
	δ,			(PRC2, KDM5A, KDM6A, KMT) --> ∅     # Degradation of histone reader and writers 
end N k0 k1 k2 d1 d0 m1 m0 p k pp kk δ v1 Kd1 v2 Kd2


# ╔═╡ 626927b6-af8f-43f5-b53d-db5ad2dab579
latexify(convert(ODESystem, Notch_model))

# ╔═╡ 7cf91626-eb96-4d42-8531-19a12f7d45e6
params(Notch_model)

# ╔═╡ a44a95ac-7cb7-4b4b-b64d-7dd43ec09629


# ╔═╡ b5078a8f-8a6c-4257-ba4a-0f3ca0fac4c0
begin
	u0 = rand(0:100,size(Notch_model.states))
	tspan = [0,1e2]
	model_p = [N, k0, k1, k2, d1, d0, m1, m0, p, k, pp, kk, δ, v1, Kd1, v2, Kd2]
	oprob = ODEProblem(Notch_model, u0, tspan, model_p)
	osol  = solve(oprob, Tsit5());
	plot(osol, vars = [:MR, :H4, :H27, :KDM5A])
end

# ╔═╡ 4655614b-7d62-44a3-8b21-8f47ef95b0b6
Notch_model

# ╔═╡ 47e7b80c-3b19-41c2-9f15-c5d216744a9a
begin
	result =[]
	for repeat = 1:10
		u0 = rand(0:100,size(Notch_model.states))
		tspan = [0,1e4]
		model_p = [N, k0, k1, k2, d1, d0, m1, m0, δ, v1, Kd1, v2, Kd2]
		oprob = ODEProblem(Notch_model, u0, tspan, p)
		osol  = solve(oprob, Tsit5());
		push!(result, osol[end][[4,6,10]])
	end 
end

# ╔═╡ 5ec364f0-f0c6-41dc-8799-013077ab7a66
# collect(round(result,2))
[round.(i,digits = 5) for i in result]

# ╔═╡ Cell order:
# ╟─347867d4-726a-4e83-977a-fd5e77512422
# ╠═d72bc2d4-b13f-450e-a24d-a056d7b48199
# ╠═dfb1f4ba-7ded-4b87-ac7d-4974349ed12d
# ╠═58c83dca-5e6c-4654-95ba-7e791fb93de6
# ╠═42fb32c8-10a6-4e75-940b-07193d4e549b
# ╠═b5302649-d9be-47e3-9a11-d39370bb678a
# ╟─d8346212-0526-4423-be21-74cde773322d
# ╠═528975db-a04c-4e2c-8aee-be2083ea6056
# ╠═825f90df-6828-440e-b7b5-07106c1ae788
# ╟─36a50cf4-2b06-4e4b-9e66-469753e4f3b4
# ╟─8732d482-987a-40a8-9da7-8ec9523931bb
# ╠═cb3bf3b5-2f4f-4307-9bf3-82022d915d3c
# ╟─e7d396f4-ecf6-4d33-80ff-2fac41990906
# ╟─626927b6-af8f-43f5-b53d-db5ad2dab579
# ╟─7cf91626-eb96-4d42-8531-19a12f7d45e6
# ╟─e456b730-ce1a-49d9-af9f-0818fd5b9d7c
# ╠═a44a95ac-7cb7-4b4b-b64d-7dd43ec09629
# ╠═b5078a8f-8a6c-4257-ba4a-0f3ca0fac4c0
# ╠═4655614b-7d62-44a3-8b21-8f47ef95b0b6
# ╠═47e7b80c-3b19-41c2-9f15-c5d216744a9a
# ╠═5ec364f0-f0c6-41dc-8799-013077ab7a66
