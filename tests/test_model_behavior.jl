using Revise
using DifferentialEquations
using Plots; pyplot()
using Catalyst
using DataFrames, CSV
using SymbolicIndexingInterface

# Include functions
include("Functions.jl")

# Load database and model
db, p_names, initi_names, parameter_set, initial_condition, rand_idx_set, _ = loading_database()
model = Import_model(; type="pulsatile")

# Test with db_idx 49
db_idx = 49
amplitude = 50.0
T_init = 1e-10
ΔT = 100
prc2 = 0.41
T_final = 1.5 * ΔT
tspan = (0.0, T_init + T_final)

# Run simulation
u0map, pmap, p, tspan, ts, cb, sol, phase = single_solve(
    model=model, 
    db_idx=db_idx, 
    freq=0.0,  # sustained
    phase=0, 
    amplitude=amplitude, 
    T_init=T_init, 
    ΔT=ΔT, 
    tspan=tspan,
    prc2=prc2
)

# Check species values at different time points
println("\nSpecies values at different times:")
println("Species order: ", [string(s) for s in species(model)])
println("\nAt t=0 (before signal):")
for (i, s) in enumerate(species(model))
    println("  $s = $(sol.u[1][i])")
end

println("\nAt t=$(T_init + 10) (signal ON):")
t_on = T_init + 10
sol_on = sol(t_on)
for (i, s) in enumerate(species(model))
    println("  $s = $(sol_on[i])")
end

println("\nAt t=$(T_init + ΔT + 10) (signal OFF):")
t_off = T_init + ΔT + 10
sol_off = sol(t_off)
for (i, s) in enumerate(species(model))
    println("  $s = $(sol_off[i])")
end

# Check the reaction behavior
println("\n\nReaction analysis:")
println("When signal is ON (A=$amplitude):")
println("  Forward rate R->NR: A * (1 + sign(cos(0))) = $amplitude * 2 = $(amplitude * 2)")
println("  Backward rate NR->R: 1.0")
println("  Net effect: R should DECREASE, NR should INCREASE")
println("  This means less free R available to bind M")
println("  So MR should DECREASE")

# Plot all species to understand the dynamics
plt = plot(sol, 
    vars=1:11,
    lw=2,
    xlabel="Time",
    ylabel="Concentration",
    title="All Species Dynamics",
    legend=:outerright,
    size=(1000, 600)
)

# Add vertical lines for signal on/off
vline!([T_init], label="Signal ON", color=:green, lw=2, ls=:dash)
vline!([T_init + ΔT], label="Signal OFF", color=:red, lw=2, ls=:dash)

savefig(plt, "test_all_species.png")
println("\nSaved plot as test_all_species.png")

# Now let's specifically check MR, H4, H27
println("\n\nFocusing on MR, H4, H27:")
println("MR index: 4")
println("H4 index: 6") 
println("H27 index: 9")

plt2 = plot(sol,
    idxs=[4, 6, 9],
    lw=2,
    xlabel="Time",
    ylabel="Concentration", 
    label=["MR" "H4" "H27"],
    title="MR, H4, H27 Dynamics"
)

# Add vertical lines for signal on/off
vline!([T_init], label="Signal ON", color=:green, lw=2, ls=:dash)
vline!([T_init + ΔT], label="Signal OFF", color=:red, lw=2, ls=:dash)

savefig(plt2, "test_MR_H4_H27.png")
println("Saved plot as test_MR_H4_H27.png")