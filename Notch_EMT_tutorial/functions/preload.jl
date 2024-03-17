using Revise
using Pkg
Pkg.activate(dirname(@__DIR__))

using DifferentialEquations
using Catalyst
using Catalyst: parameters

using ProgressMeter
using DataFrames, CSV, Random, Distributions
using StatsPlots
using Latexify, Measures, FLoops, LaTeXStrings
using LaTeXStrings


# println("Number of working processes: ", nworks())
cd(dirname(@__DIR__)); @info "Current working directory: $(pwd())"

includet(joinpath(@__DIR__, "Functions.jl"))
db, p_names, initi_names, parameter_set, initial_condition, rand_idx_set, db_idx = loading_database(; data_path="../../../../Notch_EMT_data/Notch_params_complete.csv");
function show_avail_functions()
    show(stdout, names(Main))
end
