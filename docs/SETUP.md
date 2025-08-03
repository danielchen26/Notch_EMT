# Setup Guide

This guide provides detailed instructions for setting up the Notch_EMT project environment.

## Prerequisites

- Julia 1.8 or higher
- Git
- Python (optional, for PyPlot backend)

## Installation Steps

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/Notch_EMT.git
cd Notch_EMT
```

### 2. Set Up Julia Environment

```bash
# Start Julia with the project environment
julia --project=.

# In the Julia REPL, instantiate the project
julia> using Pkg
julia> Pkg.instantiate()
```

This will install all required packages specified in `Project.toml`.

### 3. Verify Installation

Run the test scripts to verify everything is working:

```julia
# Test data access
include("tests/test_data_access.jl")

# Test paper script structure
include("tests/test_paper_structure.jl")
```

### 4. Set Up Parameter Database

The project requires a parameter database file located at `../Notch_EMT_data/Notch_params_complete.csv` (relative to the project root).

If you don't have this file:
1. Contact the project maintainers for the database
2. Or generate it using scripts in `parameter_generation/`

### 5. Optional: PyPlot Configuration

If you encounter issues with plotting:

```julia
# Rebuild PyPlot
using Pkg
Pkg.build("PyPlot")

# Set backend explicitly
using Plots
pyplot()
```

## Common Issues

### Sundials Error
If you get a Sundials library error:
```julia
using Pkg
Pkg.build("Sundials")
```

### Memory Issues
For large parameter sweeps, consider:
- Reducing the resolution of amplitude/frequency ranges
- Using pre-computed data when available
- Running on a machine with more RAM

### Missing Data Files
The following pre-computed data files speed up analysis:
- `df_592_3d.csv` - 3D parameter sweep data
- `df_49_592.csv` - Comparison data

These will be generated automatically if not found, but computation may take hours.

## Environment Variables

No special environment variables are required. The project uses relative paths for all data access.

## Development Setup

For development, also install:
```julia
using Pkg
Pkg.add("Revise")  # For interactive development
```

## Next Steps

After setup:
1. Run the tutorial: `include("Notch_EMT_tutorial/Notch_EMT_main.jl")`
2. Generate paper figures: `include("src/Notch_EMT_paper.jl")`
3. Explore example notebooks in `examples/notebooks/`