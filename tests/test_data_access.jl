# Test to verify data access for paper script
println("Testing data access for Notch_EMT_paper.jl...")
println("=" ^ 60)
println("Note: This test should be run from the project root, not from tests/")
println()

# Test 1: Check if the parameter database exists at expected location
# Adjust path based on where we're running from
if basename(pwd()) == "tests"
    data_path = "../../Notch_EMT_data/Notch_params_complete.csv"
else
    data_path = "../Notch_EMT_data/Notch_params_complete.csv"
end
if isfile(data_path)
    println("✓ Parameter database found at: $data_path")
    filesize_mb = filesize(data_path) / 1024 / 1024
    println("  Size: $(round(filesize_mb, digits=2)) MB")
else
    println("✗ Parameter database NOT found at: $data_path")
    println("  This will cause the paper script to fail!")
end

# Test 2: Check if Functions.jl can be loaded
println("\n2. Testing Functions.jl import...")
functions_path = "../src/Functions.jl"
if isfile(functions_path)
    println("✓ Functions.jl found at: $functions_path")
else
    println("✗ Functions.jl NOT found at: $functions_path")
end

# Test 3: Check figure output directories
println("\n3. Checking figure output directories...")
figure_dirs = [
    "../figures/paper",
    "../Data/regular"
]

for dir in figure_dirs
    if isdir(dir)
        println("✓ Directory exists: $dir")
    else
        println("ℹ Directory will be created: $dir")
    end
end

# Test 4: Check for pre-computed data files
println("\n4. Checking for pre-computed data files...")
precomputed_files = [
    "../df_592_3d.csv",
    "../df_49_592.csv",
    "../Data/regular/db_idx:49/pulsatile/freq :0.0:0.02:2.0_|_amplitude :0:1:300.csv"
]

for file in precomputed_files
    if isfile(file)
        println("✓ Found: $(basename(file))")
    else
        println("ℹ Not found (will be generated): $(basename(file))")
    end
end

println("\n" * "=" ^ 60)
println("Data access verification complete!")