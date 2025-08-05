# ModelingToolkit Migration Guide

This document comprehensively covers the migration from older Catalyst.jl/ModelingToolkit.jl versions to the newer API, based on issues encountered and resolved during the Notch-EMT project update.

## Overview of Changes

The main breaking changes when updating from older Catalyst/ModelingToolkit versions:
1. Parameter indexing changed from integer-based to symbolic indexing
2. MTKParameters replaced simple parameter arrays
3. SymbolicIndexingInterface is now required for parameter access in callbacks

## Key Issues and Solutions

### 1. Parameter Indexing in Callbacks

**Old Code (Integer Indexing):**
```julia
function make_cb(ts_in, index, value)
    condition(u, t, integrator) = t in ts
    function affect!(integrator)
        integrator.p[index] = value  # Integer indexing no longer works
    end
    cb = DiscreteCallback(condition, affect!)
end
```

**New Code (Symbolic Indexing):**
```julia
using SymbolicIndexingInterface  # Required import

function make_cb(ts_in, prob, model, value)
    condition(u, t, integrator) = t in ts
    # Get the parameter from the model
    A_param = parameters(model)[1]  # Get symbolic parameter
    # Get the parameter index using SymbolicIndexingInterface
    A_idx = parameter_index(prob, A_param)
    
    function affect!(integrator)
        integrator.p[A_idx] = value  # Use parameter index
    end
    cb = DiscreteCallback(condition, affect!)
end
```

### 2. Parameter Vector Construction

**Critical Issue:** The order of parameters in the model may differ from the order in your database.

**Model Parameter Order (from Catalyst):**
```julia
parameters(model) = [A, w, ϕ, k1, k2, k0, d, m, p, kk, pp, k, δ, α1]
```

**Database Parameter Order:**
```julia
database_params = [k0, k1, k2, d, m, p, k, pp, kk, δ, α1]
```

**Solution - Explicit Mapping:**
```julia
function single_solve(; model, db_idx, freq, phase, amplitude, ...)
    db_params = parameter_set[db_idx, :]
    
    # Create parameter vector in the order expected by the model
    p = [
        0.0,           # A (will be set by callback)
        freq,          # w  
        phase,         # ϕ
        db_params.k1,  # k1
        db_params.k2,  # k2
        db_params.k0,  # k0
        db_params.d,   # d
        db_params.m,   # m
        db_params.p,   # p
        db_params.kk,  # kk
        db_params.pp,  # pp
        db_params.k,   # k
        db_params.δ,   # δ
        db_params.α1   # α1
    ]
    pmap = parameters(model) .=> p
    # ... rest of function
end
```

### 3. Database Loading Issues

**Problem:** Multiple parameter database files with different formats.

**Solution:** Ensure you're loading from the correct database file that contains both parameters and initial conditions:
```julia
function loading_database(; data_path="../Notch_EMT_data/Notch_params_complete.csv")
    db = CSV.File(data_path) |> DataFrame
    
    # Parameter names are in columns 1-11
    p_names = names(db)[1:11]
    parameter_set = db[:, p_names]
    
    # Initial condition names are in columns 13-end
    initi_names = names(db)[13:end]
    initial_condition = db[:, initi_names]
    
    # ... rest of function
end
```

### 4. Signal Visualization Bug

**Problem:** Sustained signal was visualized incorrectly, showing as OFF during the ON period.

**Old Code:**
```julia
plot!(plt, [0, ts[1], ts[2], tspan[end]], [0, 0, 2 * amplitude, 0], ...)
```

**Fixed Code:**
```julia
plot!(plt, [0, ts[1], ts[1], ts[2], ts[2], tspan[end]], 
           [0, 0, 2 * amplitude, 2 * amplitude, 0, 0], ...)
```

## Debugging Strategies

### 1. Verify Model Behavior
Create a simple test to verify the model dynamics:
```julia
# Test that MR decreases when signal is ON
sol = single_solve(...)
println("MR before signal: ", sol.u[1][4])
println("MR during signal: ", sol(T_init + 10)[4])
```

### 2. Check Parameter and Species Order
```julia
println("Species order: ", species(model))
println("Parameters order: ", parameters(model))
```

### 3. Compare Against Expected Results
- Always check that biological behavior makes sense
- MR should decrease when Notch signal is ON (RBPJ sequestration)
- Use test plots to verify dynamics before generating final figures

## Best Practices for Future-Proofing

1. **Always use symbolic indexing** instead of integer indexing for parameters
2. **Document parameter order** explicitly in your code
3. **Create helper functions** for parameter mapping to centralize the logic
4. **Test model behavior** independently before complex visualizations
5. **Keep database paths configurable** to easily switch between different parameter sets

## Required Package Imports

```julia
using Catalyst
using DifferentialEquations
using SymbolicIndexingInterface  # New requirement
using DataFrames, CSV
using Plots
```

## Common Error Messages and Solutions

1. **`MethodError: no method matching setindex!(::MTKParameters{...}, ::Float64, ::Int64)`**
   - Solution: Use symbolic indexing with `parameter_index()`

2. **`ArgumentError: column name "k1" not found in the data frame`**
   - Solution: Check you're loading from the correct database file

3. **Figure shows incorrect dynamics**
   - Solution: Verify parameter mapping and signal visualization

## Testing Checklist

Before considering the migration complete:
- [ ] All callbacks use symbolic parameter indexing
- [ ] Parameter mapping is explicit and documented
- [ ] Database loading uses the correct file path
- [ ] Signal visualizations match the actual signal timing
- [ ] Generated figures match expected biological behavior
- [ ] Model dynamics are independently verified