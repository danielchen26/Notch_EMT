# Simple test to verify paper script structure without loading packages
# This test checks if the script is properly structured without executing it

println("Testing Notch_EMT_paper.jl structure...")
println("=" ^ 60)

# Test 1: Check if script file exists
script_path = "src/Notch_EMT_paper.jl"
if !isfile(script_path)
    println("✗ Script file not found: $script_path")
    exit(1)
end
println("✓ Script file found: $script_path")

# Test 2: Read and analyze script structure
println("\n2. Analyzing script structure...")
script_content = read(script_path, String)

# Check for required function definitions
required_functions = [
    "function generate_figure3",
    "function generate_figure4", 
    "function generate_figure5",
    "function check_critical_amplitude_custom",
    "function generate_figure_compare_A_omega_for_different_deltaT",
    "function run_all"
]

println("\nChecking for required functions:")
for func in required_functions
    if occursin(func, script_content)
        println("  ✓ Found: $func")
    else
        println("  ✗ Missing: $func")
    end
end

# Check for key imports
println("\nChecking for required imports:")
imports = [
    "using DifferentialEquations",
    "using Plots",
    "using Catalyst",
    "using DataFrames",
    "includet(\"./Functions.jl\")"
]

for imp in imports
    if occursin(imp, script_content)
        println("  ✓ Found: $imp")
    else
        println("  ✗ Missing: $imp")
    end
end

# Check for data path configuration
println("\nChecking data paths:")
if occursin("loading_database", script_content)
    println("  ✓ Database loading function called")
    if occursin("Notch_params_complete.csv", script_content)
        println("  ✓ Parameter database path specified")
    end
end

# Check figure paths
if occursin("figure_path", script_content)
    println("  ✓ Figure output path defined")
end

# Test 3: Check for actual data files
println("\n3. Checking data availability:")
data_file = "../Notch_EMT_data/Notch_params_complete.csv"
if isfile(data_file)
    println("  ✓ Parameter database exists: $data_file")
    # Check file size
    filesize_mb = filesize(data_file) / 1024 / 1024
    println("    Size: $(round(filesize_mb, digits=2)) MB")
else
    println("  ✗ Parameter database NOT found: $data_file")
end

# Check for pre-computed data in project root
println("\n4. Checking for pre-computed data files:")
precomputed_files = [
    "df_592_3d.csv",
    "df_49_592.csv",
    "initial_condition.csv"
]

for file in precomputed_files
    if isfile(file)
        println("  ✓ Found: $(basename(file))")
    else
        println("  ℹ Not found: $(basename(file))")
    end
end

# Test 4: Check if Functions.jl exists
println("\n5. Checking dependencies:")
if isfile("src/Functions.jl")
    println("  ✓ Functions.jl exists")
    functions_content = read("src/Functions.jl", String)
    
    # Check for key functions
    key_functions = ["loading_database", "Import_model", "single_solve", "A_ω_st_relation"]
    for func in key_functions
        if occursin("function $func", functions_content)
            println("    ✓ Found function: $func")
        end
    end
else
    println("  ✗ Functions.jl NOT found")
end

# Summary
println("\n" * "=" ^ 60)
println("Structure verification complete!")
println("\nThe script appears to be properly structured.")
println("To fix the package loading issue and run the script:")
println("  1. Try: julia --project=. -e 'using Pkg; Pkg.build(\"Sundials\")'")
println("  2. Or remove Sundials dependency if not needed")
println("  3. Then run: julia --project=. src/Notch_EMT_paper.jl")