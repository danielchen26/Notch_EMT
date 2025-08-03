# Notch_EMT Project Structure

## Overview
This document describes the reorganized structure of the Notch_EMT project as of August 2024. The reorganization was done to improve maintainability and clarity while preserving the functionality of the paper generation scripts.

## Directory Structure

```
Notch_EMT/
├── src/                        # Main source code
│   ├── Notch_EMT_paper.jl     # Main paper figure generation script ⭐
│   ├── Functions.jl           # Core simulation and analysis functions
│   └── archive/               # Old versions of main scripts
│
├── parameter_generation/       # Parameter database generation
│   ├── Notch_param_space_size.jl
│   ├── Notch_param_run_bg.jl
│   ├── Notch_param_core_space_size.jl
│   ├── Search_params.jl
│   ├── Project.toml          # Dependencies for parameter generation
│   └── Manifest.toml
│
├── tests/                     # Test scripts
│   ├── test_paper_structure.jl
│   ├── test_paper_generation.jl
│   ├── test_paper_script.jl
│   └── test_data_access.jl
│
├── Data/                      # Simulation data (git-ignored)
│   ├── parameter_databases/   # Generated parameter sets
│   ├── regular/              # Regular simulation results
│   └── Poisson/              # Poisson simulation results
│
├── figures/                   # Generated figures (git-ignored)
│   ├── paper/                # Paper figures (Fig 3, 4, 5)
│   └── 1 gene/               # Additional analysis figures
│
├── models/                    # Model definitions
│   └── legacy/               # Old model implementations
│       ├── Notch_epi.jl
│       └── Notch_epi(catalyst_version).jl
│
├── utils/                     # Utility functions
│   └── functions.jl          # Helper functions
│
├── visualization/             # Visualization scripts
│   └── Notch_parameter_space_vis.py
│
├── docs/                      # Documentation
│   ├── PROJECT_STRUCTURE.md  # This file
│   └── reorganization/       # Reorganization documentation
│
├── Notch_EMT_tutorial/       # Tutorial and examples
│   └── functions/            # Tutorial-specific functions
│
└── Notebooks/                # Jupyter notebooks
```

## Key Files

### Critical Files (DO NOT MODIFY without careful testing)
- `src/Notch_EMT_paper.jl` - Main paper figure generation script
- `src/Functions.jl` - Core functions used by paper script
- `../Notch_EMT_data/Notch_params_complete.csv` - Parameter database (external)

### Data Files
- `df_592_3d.csv` - Pre-computed 3D parameter sweep data
- `df_49_592.csv` - Pre-computed comparison data
- `initial_condition.csv` - Initial conditions for simulations

## Workflow

1. **Parameter Generation**: Scripts in `parameter_generation/` create parameter databases
2. **Simulation**: `src/Notch_EMT_paper.jl` runs simulations using parameters
3. **Analysis**: Results are saved in `Data/` and figures in `figures/`

## Important Notes

- The `Data/` and `figures/` directories are git-ignored due to large file sizes
- Parameter database is located at `../Notch_EMT_data/` (outside project root)
- Old code from `Others/Pluto_nb` and `Others/SciML_Notch` has been removed as obsolete
- A backup tag `pre-reorganization-backup` was created before reorganization

## Running the Paper Script

```bash
# From project root
julia --project=. src/Notch_EMT_paper.jl
```

## Testing

Run tests to verify the structure:
```bash
julia --project=. tests/test_paper_structure.jl
julia --project=. tests/test_data_access.jl
```