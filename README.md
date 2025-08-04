# Notch_EMT

A comprehensive Julia-based computational framework for investigating Notch signaling dynamics in Epithelial-Mesenchymal Transition (EMT). This project explores how pulsatile versus sustained Notch signaling patterns control cell fate decisions through epigenetic modifications.

## Overview

This repository contains the complete computational framework for studying:
- **Pulsatile vs Sustained Signaling**: How different temporal patterns of Notch activation affect EMT
- **Amplitude-Frequency Relationships**: The A-ω phase space that determines cell fate switching
- **Epigenetic Control**: Role of PRC2-mediated histone modifications in stabilizing cell states
- **Parameter Space Analysis**: Systematic exploration of ~1000 biologically relevant parameter sets

## Key Features

- 🧬 **Mechanistic ODE Model**: 18-species reaction network coupling Notch signaling, EMT transcription factors, and histone modifications
- 📊 **Automated Analysis**: Complete pipeline from simulation to publication-ready figures
- 🔄 **Reproducible Research**: Pre-computed data and exact figure regeneration
- ⚡ **Optimized Performance**: Efficient parameter sweeps with caching
- 📈 **Rich Visualizations**: Both Julia and Python visualization tools

## Installation

### Prerequisites
- Julia 1.10 or later
- Python 3.8+ (for PyPlot backend)
- ~4GB RAM for full analysis
- ~2GB disk space for data

### Setup
```bash
# Clone repository
git clone <repository-url>
cd Notch_EMT

# Install Julia dependencies
julia --project=.
julia> using Pkg; Pkg.instantiate()

# Download parameter database (required)
# Contact authors or see documentation for access
```

## Quick Start

### Generate All Paper Figures
```bash
julia --project=. src/Notch_EMT_paper.jl
```
This generates 8 publication-ready figures in `figures/paper/` (~10 minutes with pre-computed data).

### Run Individual Analyses
```julia
# Start Julia with project
julia --project=.

# Load the framework
include("src/Notch_EMT_paper.jl")

# Generate specific figures
generate_figure3()  # Sustained vs pulsatile dynamics
generate_figure4()  # A-ω phase diagram  
generate_figure5()  # PRC2 rate effects

# Custom analysis
include("src/Functions.jl")
model, init = Import_model()
data = loading_database()
result = single_solve(model, data, 49; freq=0.5, amplitude=100, ΔT=150)
```

### Parameter Space Visualization
```bash
# Run Python analysis (requires matplotlib, seaborn, pandas)
python visualization/parameter_analysis.py
```

## Project Structure

```
Notch_EMT/
├── src/                          # Core source code
│   ├── Functions.jl              # Model & simulation functions
│   ├── Notch_EMT_paper.jl        # Paper figure generation
│   └── archive/                  # Development history
├── parameter_generation/         # Parameter exploration tools
├── Data/                        # Simulation data (git-ignored)
│   ├── parameter_databases/     # Bistability analysis
│   ├── precomputed/            # Pre-computed results
│   └── regular/                # Simulation outputs
├── figures/                     # Generated figures (git-ignored)
├── tests/                      # Verification scripts
├── docs/                       # Documentation
│   └── COMPREHENSIVE_DOCUMENTATION.md
├── visualization/              # Python analysis tools
└── examples/notebooks/         # Jupyter explorations
```

## Model Details

The model integrates:
- **Notch Pathway**: Ligand binding, NICD release, transcriptional activation
- **EMT Core Circuit**: Snail1/2, Zeb1, miR-34, miR-200 feedback loops
- **Epigenetic Layer**: H4K20me0 (active) and H3K27me3 (repressive) marks
- **PRC2 Complex**: Dynamic regulation of histone modifications

Key innovations:
- Discrete pulse implementation for accurate pulsatile signaling
- Bistability-preserving parameter sets from biological constraints
- Systematic A-ω phase space characterization

## Documentation

- 📚 [Comprehensive Documentation](docs/COMPREHENSIVE_DOCUMENTATION.md) - Complete technical guide
- 🚀 [Setup Guide](docs/SETUP.md) - Detailed installation instructions
- 🗂️ [Project Structure](docs/PROJECT_STRUCTURE.md) - File organization
- 🤖 [AI Assistant Notes](CLAUDE.md) - Development notes and tips

## Testing

Run the test suite to verify installation:
```bash
julia --project=. tests/test_paper_script.jl
julia --project=. tests/test_data_access.jl
```

## Citation

If you use this code in your research, please cite:
```bibtex
@article{notch_emt_2024,
  title={Temporal Dynamics of Notch Signaling in Epithelial-Mesenchymal Transition},
  author={[Authors]},
  journal={[Journal]},
  year={2024},
  doi={[DOI]}
}
```

## Contributing

We welcome contributions! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

See [Contributing Guidelines](docs/CONTRIBUTING.md) for details.

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## Contact

- **Lead Developer**: [Name] ([email])
- **Principal Investigator**: [Name] ([email])
- **Issues**: Please use GitHub Issues for bug reports and feature requests
- **Parameter Database**: Contact authors for access to the full parameter database

## Acknowledgments

- Funding: [Grant information]
- Computational Resources: [HPC/Cloud provider]
- Collaborators: [Institutions]

---

*Last updated: August 2024 | Version: 2.0 (Post-reorganization)*