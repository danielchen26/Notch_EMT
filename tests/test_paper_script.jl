# Test script to verify Notch_EMT_paper.jl functionality
# This script tests each component separately

println("Starting verification of Notch_EMT_paper.jl...")
println("=" ^ 60)

# Load the necessary components without running everything
println("\n1. Loading packages and functions...")
try
    using Revise
    using DifferentialEquations
    using Plots; pyplot()
    using Catalyst
    using Catalyst: parameters
    using DataFrames, CSV, Random, Distributions
    using StatsPlots
    using ProgressMeter
    using Latexify, Measures, LaTeXStrings
    println("✓ All packages loaded successfully")
catch e
    println("✗ Error loading packages: ", e)
    exit(1)
end

# Include Functions.jl
try
    include("../src/Functions.jl")
    println("✓ Functions.jl loaded successfully")
catch e
    println("✗ Error loading Functions.jl: ", e)
    exit(1)
end

# Test database loading
println("\n2. Testing database loading...")
try
    db, p_names, initi_names, parameter_set, initial_condition, rand_idx_set, global_db_idx = loading_database(; data_path="../../Notch_EMT_data/Notch_params_complete.csv")
    println("✓ Database loaded successfully")
    println("  - Number of parameter sets: ", nrow(db))
    println("  - Parameter names: ", length(p_names))
    println("  - Initial condition variables: ", length(initi_names))
catch e
    println("✗ Error loading database: ", e)
    exit(1)
end

# Test model import
println("\n3. Testing model import...")
try
    model_pulsatile = Import_model(; type="pulsatile")
    println("✓ Model imported successfully")
    println("  - Number of species: ", length(species(model_pulsatile)))
    println("  - Number of parameters: ", length(parameters(model_pulsatile)))
catch e
    println("✗ Error importing model: ", e)
    exit(1)
end

# Define paths (adjusted for tests directory)
figure_path = joinpath(dirname(@__DIR__), "figures", "paper_test")
data_path_regular = joinpath(dirname(@__DIR__), "Data", "regular")
data_path_main = joinpath(dirname(@__DIR__), "Data")

# Create test output directory
isdir(figure_path) || mkpath(figure_path)

# Include the figure generation functions from the paper script
println("\n4. Loading figure generation functions...")
# We need to include just the functions, not run them
include_string(Main, """
# Copy of figure generation functions from Notch_EMT_paper.jl
$(read("../src/Notch_EMT_paper.jl", String)[1:end-1000])  # Exclude the run_all() call at the end
""")

println("✓ Figure generation functions loaded")

# Test each figure generation function
println("\n5. Testing figure generation functions...")
println("-" * 40)

# Test Figure 3
println("\n5.1 Testing Figure 3 generation...")
try
    fig3 = generate_figure3()
    println("✓ Figure 3 generated successfully")
    # Check if files were created
    fig3_files = ["figure3_panelA_sustained.png", "figure3_panelB_pulsatile.png", "figure3_sustained_vs_pulsatile.png"]
    for file in fig3_files
        if isfile(joinpath(figure_path, file))
            println("  ✓ Created: $file")
        else
            println("  ✗ Missing: $file")
        end
    end
catch e
    println("✗ Error generating Figure 3: ", e)
    println("Stack trace:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end

# Test Figure 4 (this might take longer due to parameter sweeps)
println("\n5.2 Testing Figure 4 generation...")
println("Note: This may take several minutes due to parameter sweeps...")
try
    # First check if pre-computed data exists
    db_idx = 49
    amplitude_range = 0:1:300
    freq_range = 0:0.02:2
    df_save_path = joinpath(dirname(@__DIR__), "Data", "regular", "db_idx:$db_idx", "pulsatile")
    filename = joinpath(df_save_path, "freq :$freq_range" * "_|_" * "amplitude :$amplitude_range.csv")
    
    if isfile(filename)
        println("  ℹ Using pre-computed data from: $(basename(dirname(filename)))/$(basename(filename))")
    else
        println("  ℹ No pre-computed data found. This will run simulations (may take 10-30 minutes)...")
    end
    
    fig4 = generate_figure4()
    println("✓ Figure 4 generated successfully")
    # Check if files were created
    fig4_files = ["figure4_panelA_boundary.png", "figure4_panelB_ST_vs_omega.png", "figure4_A_omega_relation.png"]
    for file in fig4_files
        if isfile(joinpath(figure_path, file))
            println("  ✓ Created: $file")
        else
            println("  ✗ Missing: $file")
        end
    end
catch e
    println("✗ Error generating Figure 4: ", e)
    println("Stack trace:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end

# Test Figure 5
println("\n5.3 Testing Figure 5 generation...")
try
    # Check for pre-computed data
    data_file_path = joinpath(dirname(@__DIR__), "df_592_3d.csv")
    if isfile(data_file_path)
        println("  ℹ Using pre-computed data from: df_592_3d.csv")
    else
        println("  ℹ No pre-computed data found. This will run extensive simulations (may take hours)...")
    end
    
    fig5 = generate_figure5()
    println("✓ Figure 5 generated successfully")
    if isfile(joinpath(figure_path, "figure5_A_omega_vs_PRC2.png"))
        println("  ✓ Created: figure5_A_omega_vs_PRC2.png")
    else
        println("  ✗ Missing: figure5_A_omega_vs_PRC2.png")
    end
catch e
    println("✗ Error generating Figure 5: ", e)
    println("Stack trace:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end

println("\n" * "=" * 60)
println("Verification complete!")
println("\nGenerated figures are saved in: $figure_path")
println("\nTo run the full paper script, use:")
println("  julia src/Notch_EMT_paper.jl")
println("\nTo run this test again:")
println("  cd tests && julia test_paper_script.jl")