# Notch-EMT User Guide

## Overview

The Notch-EMT project investigates how Notch signaling dynamics control epithelial-mesenchymal transition (EMT) through histone modifications. This guide covers how to use the project for simulations, analysis, and figure generation.

## Project Structure

```
Notch_EMT/
├── src/                    # Core source code
│   ├── Functions.jl        # Simulation and analysis functions
│   ├── Notch_EMT_paper.jl  # Paper figure generation
│   └── archive/            # Old versions
├── Data/                   # Simulation data (git-ignored)
│   ├── regular/            # Regular simulation results
│   └── precomputed/        # Pre-computed results
├── figures/                # Generated figures (git-ignored)
│   ├── paper/              # Paper figures
│   └── parameter_analysis/ # Parameter analysis plots
├── tests/                  # Test scripts
├── docs/                   # Documentation
└── Notch_EMT_data/         # External parameter database
```

## Key Features

### Model System
- Notch-RBPJ signaling with pulsatile/sustained inputs
- Histone modification dynamics (H4 active, H27 repressive)
- Competition between writers (KMT, PRC2) and erasers (KDM5A, KDM6A)

### Analysis Capabilities
- Amplitude-frequency-switching time relationships
- Bistability detection and classification
- Critical amplitude determination
- PRC2 rate control analysis

## Generating Paper Figures

### All Figures at Once
```julia
include("src/Notch_EMT_paper.jl")
run_all()  # Generates Figures 3, 4, and 5
```

### Individual Figures
```julia
include("src/Notch_EMT_paper.jl")
generate_figure3()  # Sustained vs pulsatile dynamics
generate_figure4()  # A-ω switching boundary
generate_figure5()  # PRC2 rate control
```

### Expected Outputs
- Figures saved to `figures/paper/`
- Figure 3: Signal type comparison (sustained vs pulsatile)
- Figure 4: Amplitude-frequency phase diagram
- Figure 5: PRC2 rate effect on switching boundary

## Running Custom Simulations

### Single Simulation
```julia
include("src/Functions.jl")

# Load database and model
db, p_names, initi_names, parameter_set, initial_condition, rand_idx_set, db_idx = loading_database()
model = Import_model(; type="pulsatile")

# Run simulation
sol = single_solve(
    model=model,
    db_idx=49,        # Parameter set ID
    freq=0.0,         # Frequency (0 = sustained)
    phase=0.0,        # Phase
    amplitude=50.0,   # Signal amplitude
    T_init=1e-10,     # Initial time
    ΔT=100,           # Signal duration
    tspan=(0.0, 200), # Total time span
    prc2=0.41         # PRC2 rate
)
```

### Parameter Sweeps
```julia
# Amplitude-frequency sweep
df = A_ω_st_relation(
    model=model,
    db_idx=49,
    amplitude_range=0:5:300,
    freq_range=0:0.1:2,
    ΔT=100
)
```

## Key Functions Reference

### Core Functions (`src/Functions.jl`)

- **`loading_database()`** - Load parameter sets and initial conditions
- **`Import_model()`** - Create Catalyst reaction network
- **`single_solve()`** - Run single simulation
- **`single_solve_plot()`** - Run simulation and plot results
- **`A_ω_st_relation()`** - Parameter sweep for A-ω analysis
- **`check_switching()`** - Detect histone state switching
- **`find_ST()`** - Calculate switching time

### Important Parameters

- **`db_idx`** - Parameter set identifier (e.g., 49, 592)
- **`amplitude`** - Notch signal amplitude
- **`freq`** - Signal frequency (0 for sustained)
- **`ΔT`** - Signal pulse duration
- **`prc2`** - PRC2 methylation rate

## Data Management

### Pre-computed Data
The project caches simulation results in `Data/regular/db_idx:X/` to avoid re-running expensive computations.

### Parameter Database
Parameters are stored in `../Notch_EMT_data/Notch_params_complete.csv` with columns for:
- Kinetic parameters (k0, k1, k2, d, m, p, k, pp, kk, δ, α1)
- Initial conditions for all species

## Troubleshooting

### Common Issues

1. **Missing data files**: Ensure `Notch_EMT_data/` folder exists at the parent directory level
2. **PyPlot issues**: Run `pyplot()` before plotting
3. **Memory issues**: Reduce parameter sweep resolution
4. **Parameter errors**: Check parameter order matches model expectations

### Performance Tips

- Use pre-computed data when available
- Reduce amplitude/frequency resolution for faster sweeps
- Run parameter sweeps in parallel when possible

## Advanced Usage

### Custom Callbacks
```julia
# Create custom callback for parameter changes
ts, cb = make_cb([T_init, T_init + ΔT], prob, model, amplitude)
sol = solve(prob, Rosenbrock23(), callback=cb, tstops=ts)
```

### Phase Control
```julia
# Disable automatic phase reset for custom phase control
sol = single_solve(..., phase_reset=false)
```

### PRC2 Rate Studies
```julia
# Sweep PRC2 rates
df = A_ω_st_relation_prc2_range(
    model=model,
    db_idx=592,
    amplitude_range=0:5:500,
    freq_range=0:0.05:1,
    prc2_range=0.1:0.1:1,
    ΔT=100
)
```

## Next Steps

- Explore different parameter sets (db_idx values)
- Modify signal patterns (frequency, amplitude, duration)
- Analyze switching dynamics under various conditions
- Extend analysis to new parameter regimes