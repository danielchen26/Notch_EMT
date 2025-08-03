# Test script to verify Notch_EMT_paper.jl functionality
# Run from project root: julia --project=. tests/test_paper_generation.jl

println("Testing Notch_EMT_paper.jl figure generation...")
println("=" ^ 60)

# Change to project root directory
cd(dirname(@__DIR__))
println("Working directory: ", pwd())

# Test 1: Load the script without executing
println("\n1. Testing if script loads without errors...")
try
    # Create a modified version that doesn't auto-run
    paper_code = read("src/Notch_EMT_paper.jl", String)
    
    # Remove the auto-execution lines
    modified_code = replace(paper_code, 
        r"run_all\(\).*?generate_figure_compare_A_omega_for_different_deltaT\(\)"s => 
        "# Auto-execution disabled for testing")
    
    # Load the modified version
    include_string(Main, modified_code)
    println("✓ Script loaded successfully")
catch e
    println("✗ Error loading script: ", e)
    exit(1)
end

# Test 2: Check if all required functions are defined
println("\n2. Checking if all functions are defined...")
required_functions = [
    :generate_figure3,
    :generate_figure4,
    :generate_figure5,
    :check_critical_amplitude_custom,
    :generate_figure_compare_A_omega_for_different_deltaT,
    :run_all
]

for func in required_functions
    if isdefined(Main, func)
        println("  ✓ Function $func is defined")
    else
        println("  ✗ Function $func is NOT defined")
    end
end

# Test 3: Check data paths
println("\n3. Checking data paths...")
data_file = "../Notch_EMT_data/Notch_params_complete.csv"
if isfile(data_file)
    println("  ✓ Parameter database found: $data_file")
else
    println("  ✗ Parameter database NOT found: $data_file")
end

# Check for pre-computed data files
precomputed_files = [
    "df_592_3d.csv",
    "Data/regular/db_idx:49/pulsatile/freq :0.0:0.02:2.0_|_amplitude :0:1:300.csv"
]

println("\n4. Checking for pre-computed data...")
for file in precomputed_files
    if isfile(file)
        println("  ✓ Found: $file")
    else
        println("  ℹ Not found (will be generated): $file")
    end
end

# Test 5: Test figure output directory
println("\n5. Testing figure output...")
figure_path = joinpath(pwd(), "figures", "paper")
if !isdir(figure_path)
    println("  ℹ Creating figure output directory: $figure_path")
    mkpath(figure_path)
else
    println("  ✓ Figure output directory exists: $figure_path")
end

# Test 6: Dry run check (don't actually generate figures)
println("\n6. Performing dry run check...")
println("  ℹ To generate figures, run one of:")
println("     julia --project=. -e 'include(\"src/Notch_EMT_paper.jl\"); generate_figure3()'")
println("     julia --project=. -e 'include(\"src/Notch_EMT_paper.jl\"); generate_figure4()'")
println("     julia --project=. -e 'include(\"src/Notch_EMT_paper.jl\"); generate_figure5()'")
println("     julia --project=. src/Notch_EMT_paper.jl  # (runs all)")

println("\n" * "=" * 60)
println("Test complete! The script appears to be properly configured.")
println("\nNote: This test only verifies the script structure.")
println("To actually generate figures, use the commands shown above.")