# Notch_EMT

A Julia-based computational framework for investigating Notch signaling dynamics and epithelial-mesenchymal transition (EMT) through mathematical modeling and systems biology approaches.

## Overview

This project explores how pulsatile Notch signaling controls cell fate decisions through epigenetic modifications. Key findings include:
- Signal amplitude and frequency determine cell fate switching
- Histone modifications (H4/H27) create bistable switches
- PRC2 rate modulates the amplitude-frequency switching boundary

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/Notch_EMT.git
cd Notch_EMT

# Activate Julia environment
julia --project=.

# Install dependencies
julia> using Pkg; Pkg.instantiate()
```

## Quick Start

### Generate Paper Figures
```julia
# Run all paper figures (Figures 3, 4, and 5)
include("src/Notch_EMT_paper.jl")
```

### Run Interactive Tutorial
```julia
include("Notch_EMT_tutorial/Notch_EMT_main.jl")
```

## Project Structure

```
Notch_EMT/
├── src/                     # Core source code
│   ├── Notch_EMT_paper.jl  # Paper figure generation
│   └── Functions.jl        # Simulation functions
├── parameter_generation/    # Parameter database scripts
├── Data/                   # Simulation results (git-ignored)
├── figures/                # Generated figures (git-ignored)
├── tests/                  # Verification scripts
├── examples/               # Example notebooks
├── docs/                   # Documentation
└── Notch_EMT_tutorial/     # Interactive tutorial
```

See [docs/PROJECT_STRUCTURE.md](docs/PROJECT_STRUCTURE.md) for detailed documentation.

## Key Features

- **Pulsatile Signal Modeling**: Time-varying Notch input signals
- **Epigenetic Dynamics**: Histone modification state transitions
- **Parameter Sweeps**: Systematic amplitude-frequency analysis
- **Bistable Switches**: H4 (active) vs H27 (repressive) states

## Dependencies

- Julia 1.8+
- Key packages: DifferentialEquations, Catalyst, Plots, DataFrames
- Full list in Project.toml

## Citation

If you use this code, please cite:
[Paper citation to be added]

## License

[License information to be added]