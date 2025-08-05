# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Julia-based scientific computing project investigating Notch signaling dynamics and epithelial-mesenchymal transition (EMT) using differential equation models and systems biology approaches. The project focuses on how signal properties (amplitude and frequency) control cell fate switching through histone modifications.

## Key Architecture

### Core Model System
The project implements a reaction network model defined in Catalyst.jl that captures:
- Notch-RBPJ signaling dynamics with time-varying inputs (sustained vs pulsatile)
- Histone modification dynamics (H4 active marks, H27 repressive marks)
- Competition between histone writers (KMT, PRC2) and erasers (KDM5A, KDM6A)

### Data Flow
1. **Parameter Database**: Pre-computed parameter sets stored in `Notch_EMT_data/Notch_params_complete.csv`
2. **Simulation Engine**: Uses DifferentialEquations.jl with callbacks for signal control
3. **Analysis Pipeline**: Systematic parameter sweeps for amplitude-frequency-switching time relationships
4. **Visualization**: Publication-quality figures using Plots.jl with PyPlot backend

## Development Commands

### Environment Setup
```bash
# Activate Julia environment
julia --project=.

# In Julia REPL, instantiate packages
julia> using Pkg; Pkg.instantiate()
```

### Running Simulations

#### Generate Paper Figures
```julia
# Generate all figures from the paper
include("src/Notch_EMT_paper.jl")
run_all()  # Generates Figures 3, 4, and 5

# Generate individual figures
generate_figure3()  # Sustained vs pulsatile dynamics
generate_figure4()  # A-ω switching boundary
generate_figure5()  # PRC2 rate control
```

#### Tutorial/Interactive Analysis
```julia
# Run the tutorial notebook
include("Notch_EMT_tutorial/Notch_EMT_main.jl")
```

### Key Functions in `src/Functions.jl`

- `loading_database()`: Loads parameter sets and initial conditions
- `Import_model()`: Creates the Catalyst reaction network model
- `single_solve()`: Performs single simulation with specified parameters
- `A_ω_st_relation()`: Systematic parameter sweep for amplitude-frequency analysis
- `check_switching()`: Determines if histone state switching occurred
- `find_ST()`: Calculates precise switching time

## Project Structure

- `src/`: Main source code
  - `Functions.jl`: Core simulation and analysis functions (updated for new ModelingToolkit API)
  - `Notch_EMT_paper.jl`: Script to generate paper figures
  - `archive/`: Old versions of main scripts
- `tests/`: Test scripts for verification and debugging
  - Model behavior tests, parameter ordering tests, database loading tests
  - See `tests/README.md` for details
- `docs/`: Project documentation
  - `README.md`: Documentation index
  - `SETUP.md`: Installation and setup guide
  - `USER_GUIDE.md`: Comprehensive usage guide
  - `ModelingToolkit_Migration_Guide.md`: API migration guide
- `Data/`: Simulation results organized by parameter indices (git-ignored)
  - `regular/`: Regular simulation results by db_idx
  - `precomputed/`: Pre-computed analysis results
- `figures/`: Generated plots and visualizations (git-ignored)
  - `paper/`: Paper figures (Figures 3, 4, 5)
  - `parameter_analysis/`: Parameter analysis plots
- `parameter_generation/`: Scripts for generating parameter databases
- `Notch_EMT_tutorial/`: Interactive tutorial with example analyses
- `examples/`: Example analyses and notebooks
- `utils/`: Utility functions
- `visualization/`: Visualization scripts
- `models/legacy/`: Old model implementations

## Important Technical Notes

1. **Time Parameters**: 
   - `T_init`: Must be > 0 (typically 1e-10) for numerical stability
   - `ΔT`: Signal duration (default 100)
   - Phase reset: Automatic phase adjustment for pulsatile signals to ensure consistent initial conditions

2. **Plotting Backend**: Uses PyPlot for consistency. Call `pyplot()` before generating figures.

3. **Data Caching**: The scripts check for pre-computed data files before running expensive simulations. Results are saved in `Data/regular/db_idx:X/` directories.

4. **Parameter Identifiers**: 
   - `db_idx`: Indexes specific parameter sets from the database
   - Common indices: 49 (used in Figures 3,4), 592 (used in Figure 5)

5. **ModelingToolkit Compatibility**: 
   - Uses SymbolicIndexingInterface for parameter access in callbacks
   - Parameter indexing is symbolic, not integer-based
   - Fixed signal visualization bug in `single_solve_plot()`
   - See `docs/ModelingToolkit_Migration_Guide.md` for details

6. **Database Path**: 
   - Parameter database location: `../Notch_EMT_data/Notch_params_complete.csv`
   - Contains both parameters (k0, k1, k2, etc.) and initial conditions
   - `loading_database()` uses default path to this file

## Debugging Tips

- If plots don't display correctly, ensure PyPlot backend is active: `pyplot()`
- For performance issues, check if pre-computed data exists in `Data/` directory
- Phase reset issues: Set `phase_reset=false` if custom phase control needed
- Memory issues with large parameter sweeps: Reduce resolution in `amplitude_range` or `freq_range`
- Parameter mismatch errors: Check parameter order matches between database and model
- Signal visualization issues: Verify signal timing in plots matches actual callback timing
- Database loading errors: Ensure `Notch_EMT_data/` folder exists at parent directory level

## Recent Updates (2024)

- Migrated to newer Catalyst.jl/ModelingToolkit.jl API
- Fixed parameter indexing to use symbolic approach
- Corrected signal visualization timing bug
- Consolidated documentation and removed redundancy
- Added comprehensive test suite in `tests/` folder