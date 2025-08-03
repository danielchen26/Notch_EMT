# Codebase Reorganization Plan for Notch_EMT Project

## Executive Summary

This document outlines a comprehensive plan to reorganize the Notch_EMT codebase for improved maintainability, clarity, and collaboration. The reorganization addresses current issues including code duplication, scattered data files, and unclear directory structure.

## Current State Analysis

### Issues Identified

1. **Code Duplication**
   - `Functions.jl` exists in both `src/` and `Notch_EMT_tutorial/functions/`
   - Potential inconsistencies between versions
   - Unclear which version is authoritative

2. **Scattered Files**
   - Data files (CSV, PNG) mixed in root directory
   - No clear separation between input data, processed data, and outputs
   - Figures not organized by purpose or origin

3. **Unclear Organization**
   - "Others" folder contains mixed content (experiments, examples, utilities)
   - No clear distinction between core code and experimental code
   - External examples mixed with project code

4. **Legacy Code**
   - `src/archive/` contains old versions (main.jl through main5_*.jl)
   - Unclear which code is actively used vs. historical

5. **File System Issues**
   - Many `._*` prefixed files (macOS metadata) throughout
   - These files clutter the repository and serve no purpose

6. **Mixed Languages**
   - Python visualization code mixed with Julia simulation code
   - No clear separation of language-specific code

### Active vs. Legacy Code Assessment

**Currently Active:**
- `src/Notch_EMT_paper.jl` - Primary script for paper figure generation
- `src/Functions.jl` - Core function library (actively modified)
- `src/main4_postive_discrete_pulse.jl` - Appears to be in use

**Tutorial/Documentation:**
- `Notch_EMT_tutorial/` - Interactive tutorial materials

**Experimental/Exploratory:**
- `Others/SciML_Notch/` - SciML experiments
- `Others/Pluto_nb/` - Pluto notebooks for exploration

**Legacy/Archive:**
- `src/archive/main*.jl` - Old versions of main analysis scripts
- Files with `._` prefix - macOS metadata files

## Proposed New Structure

```
Notch_EMT/
│
├── src/                          # Core source code
│   ├── core/                     # Core functionality
│   │   ├── Functions.jl          # Main function library (consolidated)
│   │   ├── Models.jl             # Model definitions (extracted from Functions.jl)
│   │   ├── Parameters.jl         # Parameter handling and structures
│   │   ├── Callbacks.jl          # Callback functions for simulations
│   │   └── Utils.jl              # Utility functions
│   │
│   ├── analysis/                 # Analysis scripts
│   │   ├── Notch_EMT_paper.jl   # Paper figure generation (current version)
│   │   ├── switching_analysis.jl # Switching time analysis functions
│   │   ├── parameter_sweep.jl    # Parameter sweep functions
│   │   └── boundary_analysis.jl  # A-ω boundary analysis
│   │
│   └── visualization/            # Visualization code
│       ├── plotting.jl           # Julia plotting functions
│       ├── animations.jl         # Animation generation functions
│       └── python/               # Python visualization
│           └── parameter_space_vis.py
│
├── notebooks/                    # Interactive notebooks
│   ├── tutorial/                 # Tutorial materials
│   │   ├── Notch_EMT_tutorial.jl # Main tutorial
│   │   ├── functions/            # Tutorial-specific functions
│   │   └── README.md             # Tutorial documentation
│   │
│   └── pluto/                    # Pluto notebooks
│       ├── histone_model.jl
│       ├── histone_model_rd.jl
│       └── README.md
│
├── experiments/                  # Experimental/exploratory code
│   ├── duffing_oscillator/       # Duffing oscillator experiments
│   │   ├── A_w_curve.jl
│   │   └── Duffing_oscillator.jl
│   │
│   ├── prc2_competition/         # PRC2 competition analysis
│   │   ├── PRC2_competition.jl
│   │   └── Functions.jl
│   │
│   ├── parameter_generation/     # Parameter space exploration
│   │   ├── Notch_param_core_space_size.jl
│   │   ├── Notch_param_space_size.jl
│   │   └── Search_params.jl
│   │
│   └── sciml_explorations/       # SciML experiments
│       └── notch_EMT_dd.jl
│
├── data/                         # All data files
│   ├── raw/                      # Raw input data
│   │   └── parameters/           # Parameter databases
│   │       └── Notch_params_complete.csv
│   │
│   ├── processed/                # Processed data
│   │   ├── switching_analysis/   # Switching analysis results
│   │   │   ├── df_49_592.csv
│   │   │   └── df_592_3d.csv
│   │   │
│   │   └── initial_conditions/   # Initial condition files
│   │       └── initial_condition.csv
│   │
│   └── outputs/                  # Simulation outputs
│       └── regular/              # Regular simulation results
│           └── db_idx_*/         # Organized by parameter set ID
│
├── figures/                      # All generated figures
│   ├── paper/                    # Publication-ready figures
│   │   ├── figure3_sustained_vs_pulsatile.png
│   │   ├── figure4_A_omega_relation.png
│   │   └── figure5_A_omega_vs_PRC2.png
│   │
│   ├── analysis/                 # Analysis figures
│   │   ├── switching_time/
│   │   ├── parameter_sweeps/
│   │   └── boundaries/
│   │
│   └── exploratory/              # Exploratory plots
│       └── working/              # Temporary working figures
│
├── scripts/                      # Standalone executable scripts
│   ├── run_paper_figures.jl     # Generate all paper figures
│   ├── run_full_analysis.jl     # Run complete analysis pipeline
│   ├── process_data.jl          # Data processing utilities
│   └── clean_repository.jl      # Repository cleaning script
│
├── tests/                        # Test files
│   ├── test_models.jl            # Model definition tests
│   ├── test_functions.jl         # Function library tests
│   ├── test_analysis.jl          # Analysis function tests
│   └── test_data_integrity.jl   # Data integrity checks
│
├── docs/                         # Documentation
│   ├── CLAUDE.md                 # AI assistant guidance (existing)
│   ├── API.md                    # Function API documentation
│   ├── MODEL_DESCRIPTION.md      # Model equations and parameters
│   ├── DATA_DICTIONARY.md        # Description of all data files
│   └── tutorials/                # Tutorial documentation
│       └── getting_started.md
│
├── legacy/                       # Legacy code (preserved but isolated)
│   ├── archive/                  # Old main files from src/archive/
│   ├── sciml_examples/           # External SciML examples
│   └── README.md                 # Explanation of legacy code
│
├── Project.toml                  # Julia project file
├── Manifest.toml                 # Julia manifest
├── pyproject.toml                # Python project file
├── poetry.lock                   # Python dependencies
├── README.md                     # Project readme (updated)
├── REORGANIZATION_PLAN.md        # This document
├── CHANGELOG.md                  # Track changes
└── .gitignore                    # Comprehensive gitignore
```

## Implementation Plan

### Phase 1: Preparation and Planning (Current)
- [x] Create `paper-revision` branch
- [x] Document reorganization plan (this document)
- [ ] Review plan with stakeholders
- [ ] Create migration scripts

### Phase 2: Core Structure Setup (Days 1-2)
- [ ] Create new directory structure
- [ ] Set up comprehensive `.gitignore`
- [ ] Create placeholder README files in each directory

### Phase 3: Code Consolidation (Days 3-5)
- [ ] Merge and consolidate `Functions.jl` files
  - Compare versions for differences
  - Merge unique functions
  - Resolve conflicts
  - Document deprecated functions
- [ ] Extract model definitions into `Models.jl`
- [ ] Extract parameter structures into `Parameters.jl`
- [ ] Create `Utils.jl` for general utilities

### Phase 4: File Migration (Days 6-7)
- [ ] Move active source files to new locations
- [ ] Relocate data files to appropriate directories
- [ ] Organize figures by category
- [ ] Move notebooks to dedicated directory
- [ ] Archive legacy code

### Phase 5: Code Updates (Days 8-10)
- [ ] Update all import paths in Julia files
- [ ] Update data loading paths
- [ ] Update figure saving paths
- [ ] Test all major functions
- [ ] Update documentation

### Phase 6: Cleanup (Day 11)
- [ ] Remove all `._*` files
- [ ] Remove empty directories
- [ ] Clean up root directory
- [ ] Update README.md

### Phase 7: Validation (Day 12)
- [ ] Run all paper figure generation scripts
- [ ] Verify tutorial still works
- [ ] Run any existing tests
- [ ] Document any breaking changes

## Migration Script Outline

```julia
# migration_script.jl
# This script will help automate the reorganization

# 1. Create new directory structure
function create_directories()
    # Implementation here
end

# 2. File mapping dictionary
file_mappings = Dict(
    "src/Functions.jl" => "src/core/Functions.jl",
    "src/Notch_EMT_paper.jl" => "src/analysis/Notch_EMT_paper.jl",
    # ... more mappings
)

# 3. Move files according to mapping
function migrate_files(mappings)
    # Implementation here
end

# 4. Update import statements
function update_imports(file_path)
    # Implementation here
end

# 5. Clean up macOS metadata files
function clean_metadata_files()
    # Implementation here
end
```

## Risk Assessment and Mitigation

### Risks:
1. **Breaking existing workflows**
   - Mitigation: Create migration guide, update all paths systematically
   
2. **Losing track of file history**
   - Mitigation: Use `git mv` for all file moves to preserve history
   
3. **Import path errors**
   - Mitigation: Systematic search and replace, comprehensive testing
   
4. **Data file references**
   - Mitigation: Create path configuration file

### Benefits:
1. **Improved maintainability**: Clear structure makes code easier to navigate
2. **Better collaboration**: Contributors can easily find relevant code
3. **Reduced duplication**: Single source of truth for all functions
4. **Cleaner repository**: Organized data and output files
5. **Easier testing**: Dedicated test directory encourages test writing
6. **Better documentation**: Centralized documentation location

## Decision Points

Before proceeding, please consider:

1. **Timing**: Is now the right time for this reorganization?
2. **Scope**: Should we do the full reorganization or start with high-priority items?
3. **Automation**: Should we create migration scripts or do it manually?
4. **Testing**: Do we need to create tests before reorganizing?
5. **Documentation**: Should we document the current API before moving files?

## Next Steps

1. Review this plan and provide feedback
2. Decide on implementation timeline
3. Determine if any modifications are needed
4. Choose between automated vs. manual migration
5. Set up any necessary tooling or scripts

## Questions for Consideration

1. Are there any critical workflows that must not be disrupted?
2. Are there any files or directories that should not be moved?
3. Do you have preferences for the naming conventions?
4. Should we preserve the "Others" directory in any form?
5. How important is it to preserve the exact git history for moved files?

---

This reorganization plan is designed to transform the Notch_EMT repository into a well-organized, maintainable codebase suitable for both research and potential public release. The structure follows best practices for scientific computing projects and will make collaboration and future development much easier.