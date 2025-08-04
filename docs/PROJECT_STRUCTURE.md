# Project Structure

This document provides a detailed overview of the Notch_EMT project organization after the August 2024 reorganization.

**Last Updated**: August 2024  
**Status**: Current (Post-reorganization)

## Directory Layout

```
Notch_EMT/
├── src/                              # Core source code
│   ├── Functions.jl                  # Core simulation and analysis functions
│   ├── Notch_EMT_paper.jl           # Main paper figure generation script ⭐
│   ├── main4_postive_discrete_pulse.jl  # Additional pulse analysis
│   └── archive/                      # Historical development versions
│       ├── main.jl                   # Original development script
│       ├── main2.jl                  # A-ω collection version
│       ├── main3.jl                  # Phase animation version
│       ├── main4_*.jl                # Various main4 iterations
│       └── main5_*.jl                # Modularized version attempts
│
├── parameter_generation/             # Parameter space exploration
│   ├── Notch_param_space_size.jl    # Parameter space sizing
│   ├── Notch_param_core_space_size.jl  # Core parameter analysis
│   ├── Notch_param_run_bg.jl        # Background parameter runs
│   ├── Search_params.jl             # Parameter search utilities
│   ├── Project.toml                 # Separate environment
│   └── Manifest.toml                # Locked dependencies
│
├── Data/                            # Simulation data (git-ignored)
│   ├── parameter_databases/         # Parameter sets and analysis
│   │   ├── Notch_core_bistability_db.csv     # 868,649 parameters
│   │   ├── Notch_parameter_bistability_r|nr_statistics.csv
│   │   ├── db.csv                   # Selected parameter sets
│   │   └── db_complete.csv          # Complete database
│   │
│   ├── precomputed/                # Pre-computed results
│   │   ├── df_592_3d.csv           # PRC2 sweep (Fig 5)
│   │   ├── df_49_592.csv           # Combined analysis
│   │   ├── df_43_52.csv            # Additional analysis
│   │   └── initial_condition.csv   # Standard initial states
│   │
│   ├── regular/                    # Organized by parameter ID
│   │   └── db_idx:*/               # Results for each parameter
│   │       └── *.csv               # A-ω sweep results
│   │
│   └── Poisson/                    # Stochastic simulations
│       └── *.csv                   # Poisson process results
│
├── figures/                         # Generated figures (git-ignored)
│   ├── paper/                      # Publication figures
│   │   ├── figure3_panelA_sustained.png     # Sustained signaling
│   │   ├── figure3_panelB_pulsatile.png     # Pulsatile examples
│   │   ├── figure3_sustained_vs_pulsatile.png # Combined
│   │   ├── figure4_panelA_boundary.png      # A-ω heatmap
│   │   ├── figure4_panelB_ST_vs_omega.png   # Switching time
│   │   ├── figure4_A_omega_relation.png     # Combined
│   │   ├── figure5_A_omega_vs_PRC2.png      # PRC2 effects
│   │   └── figure_compare_A_omega_deltaT_db49.png
│   │
│   └── parameter_analysis/         # Parameter visualizations
│       ├── stability_distribution.png
│       ├── parameter_distributions.png
│       ├── parameter_violins.png
│       ├── correlation_heatmap.png
│       └── scatter_matrix.png
│
├── tests/                          # Verification scripts
│   ├── test_paper_script.jl       # Comprehensive tests
│   ├── test_paper_generation.jl   # Figure generation checks
│   ├── test_paper_structure.jl    # Structure verification
│   └── test_data_access.jl        # Data availability tests
│
├── docs/                           # Documentation
│   ├── COMPREHENSIVE_DOCUMENTATION.md # Complete guide
│   ├── PROJECT_STRUCTURE.md         # This file
│   ├── SETUP.md                     # Installation guide
│   └── reorganization/              # Reorganization history
│       ├── REORGANIZATION_PLAN.md
│       └── REORGANIZATION_SUMMARY.md
│
├── examples/notebooks/             # Jupyter explorations
│   ├── Models.ipynb               # Model exploration
│   ├── Models_stability.ipynb     # Stability analysis
│   └── A_w_curve.ipynb           # A-ω relationships
│
├── models/legacy/                  # Previous model versions
│   ├── Notch_epi.jl               # Original model
│   ├── Notch_epi(catalyst_version).jl # Catalyst port
│   └── Notch_epi_rd.ipynb        # Reaction-diffusion
│
├── utils/                         # Utility functions
│   └── functions.jl              # Helper utilities
│
├── visualization/                 # Python analysis tools
│   └── parameter_analysis.py     # Parameter visualization
│
├── Project.toml                  # Julia dependencies
├── Manifest.toml                 # Locked versions
├── README.md                     # Project overview
├── CLAUDE.md                     # AI assistant notes
└── .gitignore                    # Excludes Data/, figures/, .venv/
```

## Key Files and Their Roles

### Core Functionality
1. **`src/Functions.jl`** (2000+ lines)
   - `Import_model()`: Catalyst model definition
   - `single_solve()`: ODE simulation engine
   - `A_ω_st_relation()`: Parameter sweep analysis
   - `loading_database()`: Parameter loading
   - Plotting and analysis utilities

2. **`src/Notch_EMT_paper.jl`** (600+ lines)
   - `generate_figure3()`: Sustained vs pulsatile
   - `generate_figure4()`: A-ω phase diagram
   - `generate_figure5()`: PRC2 rate effects
   - `check_critical_amplitude_custom()`: Boundary detection
   - Complete paper reproduction pipeline

### Data Files
1. **Parameter Database** (External)
   - Location: `../Notch_EMT_data/Notch_params_complete.csv`
   - Contains: ~1000 biologically constrained parameter sets
   - Required for all simulations

2. **Pre-computed Results**
   - `df_592_3d.csv`: 592 parameters × 10 PRC2 rates × A-ω grid
   - `df_49_592.csv`: Comparative analysis results
   - Enable quick figure regeneration

### Configuration
- **Project.toml**: Package requirements
  - DifferentialEquations v7.13+
  - Catalyst v14.3+
  - Plots (PyPlot backend)
  - DataFrames, CSV, etc.

## Data Organization

### Parameter ID System
- Each parameter set has unique ID (e.g., db_idx:49)
- Results stored in corresponding directories
- Enables parallel parameter exploration

### File Naming Convention
```
freq :0.0:0.02:2.0_|_amplitude :0:1:300.csv
```
- Describes sweep parameters
- Preserves analysis conditions
- Enables result reconstruction

## Workflow Paths

### 1. Paper Figure Generation
```
Parameters → Functions.jl → Simulations → Data/regular/ → Plots → figures/paper/
```

### 2. Parameter Exploration
```
parameter_generation/ → Bistability Test → parameter_databases/ → Analysis
```

### 3. Custom Analysis
```
Import_model() → Custom Parameters → single_solve() → Analysis → Visualization
```

## Important Notes

### Git-Ignored Directories
- `Data/`: Large simulation outputs (~1GB+)
- `figures/`: Generated plots
- `.venv/`: Python virtual environment

### External Dependencies
- Parameter database must be at `../Notch_EMT_data/`
- Python required for PyPlot backend
- ~4GB RAM for full analysis

### Removed Components
- `Notch_EMT_tutorial/`: Redundant with paper code
- `Others/Pluto_nb/`: Outdated notebooks
- `Others/SciML_Notch/`: Experimental neural ODE
- Poetry files: Unused Python configuration

## Testing and Verification

### Run All Tests
```bash
julia --project=. tests/test_paper_script.jl
```

### Check Structure
```bash
julia --project=. tests/test_paper_structure.jl
```

### Verify Data Access
```bash
julia --project=. tests/test_data_access.jl
```

## Maintenance Guidelines

1. **Before Modifying Core Files**:
   - Run all tests
   - Create feature branch
   - Document changes

2. **Adding New Analysis**:
   - Implement in Functions.jl
   - Add test coverage
   - Update documentation

3. **Data Management**:
   - Use consistent naming
   - Document parameters
   - Keep backups of pre-computed data

---

*This structure represents the cleaned and organized state after removing obsolete components and consolidating functionality.*