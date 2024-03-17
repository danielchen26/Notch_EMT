@info "Loading bump model "
signal_type = "bump"
model_bump = Import_model(; type=signal_type);
equations(model_bump); @nl


@info "Loading pulsatile model "
signal_type = "pulsatile"
model_pulsatile = Import_model(; type=signal_type);
equations(model_pulsatile); @nl
rn_latex, ode_latex = ode_latex_gen(model_pulsatile);

# ====== print model vars and params =======
# @info species(model)
# @info parameters(model)
# @info speciesmap(model)
# @info paramsmap(model)

# ====== print CRN and ODE equations =====
# latexify(model_pulsatile)
# # Graph(model_pulsatile)
# ODE_equations = convert(ODESystem, model_pulsatile)
# @show latexify(ODE_equations)

# println("The CRN for pulsatile model is: $rn_latex")