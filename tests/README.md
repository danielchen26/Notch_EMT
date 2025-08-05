# Tests Directory

This directory contains test scripts for the Notch-EMT project. These tests were created during debugging and migration to newer Catalyst.jl/ModelingToolkit.jl versions.

## Test Scripts

### Core Functionality Tests

- **`test_model_behavior.jl`** - Tests the basic model dynamics to verify that MR decreases when signal is ON (critical for biological correctness)
- **`test_species_order.jl`** - Verifies the order of species and parameters in the model matches the database
- **`test_db_loading.jl`** - Tests loading from different database files to ensure correct parameter file is used

### Paper Generation Tests

- **`test_paper_script.jl`** - Tests the main paper script execution
- **`test_paper_generation.jl`** - Tests individual figure generation functions
- **`test_paper_structure.jl`** - Verifies the structure and organization of paper generation code
- **`test_data_access.jl`** - Tests data loading and access patterns

## Running Tests

To run any test script:
```julia
julia --project=. tests/test_script_name.jl
```

## Notes

These test scripts were instrumental in identifying and fixing:
1. Parameter indexing issues with new ModelingToolkit API
2. Database file path problems
3. Signal visualization bugs
4. Parameter ordering mismatches

They can be used as reference for debugging similar issues in the future.