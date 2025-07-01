# Analysis: PRC2 Boundary Crossing Limitation for Sequential Gene Activation

## Current Limitation Identified

### **Problem Statement**
The current chemical reaction network model has a **fundamental limitation**: it only supports **monotonically increasing boundary curves** when PRC2 rate changes, but **cannot model scenarios where PRC2 rate changes enable boundary crossing that leads to sequential activation of two genes in a specific order**.

### **Current Model Architecture Issues**

#### 1. **Single PRC2 Rate Per Simulation**
```julia
function remake_prob(model, u0map, tspan, p; prc2=0.4, mute_parameter_disp=false)
    # Current implementation: ONLY allows single PRC2 value
    p[5] = prc2  # Parameter index 5 is PRC2 rate - FIXED throughout simulation
    pmap = parameters(model) .=> p
    prob_new = ODEProblem(model, u0map, tspan, pmap)
end
```

**Issue**: PRC2 rate is **static** during each simulation. The model cannot handle:
- **Dynamic PRC2 changes during simulation**
- **Gene-specific PRC2 rates**
- **Temporal PRC2 rate modulation**

#### 2. **Independent Gene Treatment**
```julia
function Two_Genes_TS_by_Prc2(; model=model, id1=592, id2=49, ...)
    # Gene 1 simulation
    sol_gene1, ts1, tspan1 = remake_solve_prob_by_ID(; model=model, db_idx=id1, prc2=prc2)
    
    # Gene 2 simulation - SAME PRC2 rate
    sol_gene2, ts2, tspan2 = remake_solve_prob_by_ID(; model=model, db_idx=id2, prc2=prc2)
end
```

**Issue**: Both genes use the **same PRC2 rate** and are simulated **independently**. This prevents:
- **Sequential activation dependencies**
- **Gene 1 → Gene 2 signaling cascades**
- **PRC2-mediated cross-gene regulation**

#### 3. **Monotonic Boundary Constraint**
Current A-ω boundary generation:
```julia
function A_ω_st_relation_prc2_range(; prc2_range="NA", ...)
    for prc2 in prc2_range
        # Each PRC2 value generates independent boundary
        # No mechanism for boundary CROSSING or NON-MONOTONIC behavior
    end
end
```

**Issue**: Each PRC2 rate generates a **separate boundary curve**. The model cannot represent:
- **Boundary intersections**
- **Non-monotonic relationships**
- **Critical PRC2 transitions that flip gene activation order**

### **Missing Biological Phenomena**

#### 1. **PRC2-Mediated Gene Ordering**
**Real Biology**: PRC2 levels can determine which gene activates first in a cascade:
- **Low PRC2**: Gene A activates before Gene B
- **High PRC2**: Gene B activates before Gene A  
- **Critical PRC2**: Boundary crossing point where ordering switches

#### 2. **Sequential Activation Cascades**
**Real Biology**: Gene 1 activation → produces signals → influences Gene 2 activation
- **Feed-forward loops**
- **Temporal delays between gene activations**
- **PRC2-dependent timing control**

#### 3. **Competitive Chromatin Dynamics**
**Real Biology**: Limited PRC2 pools create competition between genes:
- **PRC2 titration effects**
- **Chromatin domain competition**
- **Context-dependent PRC2 recruitment**

## Proposed Solutions

### **Solution 1: Dynamic PRC2 Model**
Implement time-dependent PRC2 that can change during simulation:

```julia
# NEW: Time-dependent PRC2 function
function dynamic_prc2(t, base_prc2, modulation_params)
    # Allow PRC2 to change over time
    if t < modulation_params.switch_time
        return base_prc2
    else
        return modulation_params.new_prc2
    end
end

# Modified reaction network with dynamic PRC2
model_dynamic_prc2 = @reaction_network begin
    # ... existing reactions ...
    # NEW: PRC2 rate changes based on gene states or time
    dynamic_prc2(t, base_prc2, modulation), H0 + PRC2 --> H27 + PRC2
end
```

### **Solution 2: Coupled Two-Gene Model**
Create a single model with both genes and shared PRC2 pool:

```julia
model_coupled_genes = @reaction_network begin
    # Gene 1 reactions
    (A1 * signal1(t), 1.0), R1 ↔ NR1
    (k1, k2), M1 + R1 ↔ MR1
    k0, MR1 --> MR1 + KDM5A
    
    # Gene 2 reactions  
    (A2 * signal2(t), 1.0), R2 ↔ NR2
    (k1, k2), M2 + R2 ↔ MR2
    k0, MR2 --> MR2 + KDM5A
    
    # SHARED PRC2 pool - creates competition
    m, H0_gene1 + PRC2 --> H27_gene1 + PRC2
    m, H0_gene2 + PRC2 --> H27_gene2 + PRC2
    
    # Cross-gene signaling (NEW)
    signal_strength, H4_gene1 --> H4_gene1 + Gene1_Signal
    activation_rate, Gene1_Signal + R2 --> NR2 + Gene1_Signal
    
    # PRC2 production dependent on gene states (NEW)
    prc2_prod_gene1, H27_gene1 --> H27_gene1 + PRC2
    prc2_prod_gene2, H27_gene2 --> H27_gene2 + PRC2
end
```

### **Solution 3: Non-Monotonic Boundary Analysis**
Implement boundary crossing detection:

```julia
function analyze_boundary_crossing(prc2_range, amplitude_range, freq_range)
    boundaries = []
    
    for prc2 in prc2_range
        boundary = extract_min_amp(A_ω_st_relation(prc2=prc2))
        push!(boundaries, (prc2=prc2, boundary=boundary))
    end
    
    # NEW: Detect crossing points
    crossing_points = find_boundary_crossings(boundaries)
    gene_order_switches = analyze_activation_order_changes(crossing_points)
    
    return boundaries, crossing_points, gene_order_switches
end

function find_boundary_crossings(boundaries)
    crossings = []
    for i in 1:length(boundaries)-1
        intersection = find_curve_intersection(boundaries[i], boundaries[i+1])
        if !isnothing(intersection)
            push!(crossings, intersection)
        end
    end
    return crossings
end
```

### **Solution 4: Sequential Activation Framework**
Implement explicit gene ordering analysis:

```julia
function analyze_sequential_activation(prc2_val, gene1_params, gene2_params)
    # Simulate coupled system
    sol_coupled = solve_coupled_genes(prc2=prc2_val, gene1=gene1_params, gene2=gene2_params)
    
    # Extract activation times
    t_activation_gene1 = find_activation_time(sol_coupled, gene=1)
    t_activation_gene2 = find_activation_time(sol_coupled, gene=2)
    
    # Determine activation order
    if t_activation_gene1 < t_activation_gene2
        return "Gene1_first", (t_activation_gene1, t_activation_gene2)
    else
        return "Gene2_first", (t_activation_gene1, t_activation_gene2)
    end
end

function find_prc2_order_switch_point(gene1_params, gene2_params, prc2_range)
    orders = []
    for prc2 in prc2_range
        order, times = analyze_sequential_activation(prc2, gene1_params, gene2_params)
        push!(orders, (prc2=prc2, order=order, times=times))
    end
    
    # Find where order switches
    switch_point = find_order_transition(orders)
    return switch_point
end
```

## Implementation Priority

### **Immediate Implementation (High Priority)**
1. **Coupled Two-Gene Model**: Single simulation with shared PRC2 pool
2. **Cross-Gene Signaling**: Gene 1 → Gene 2 activation pathways
3. **Sequential Activation Detection**: Measure and compare activation times

### **Medium-Term Implementation**
1. **Dynamic PRC2 Modulation**: Time-dependent PRC2 changes
2. **Boundary Crossing Analysis**: Non-monotonic boundary relationships
3. **Gene Competition Effects**: Limited PRC2 pool dynamics

### **Advanced Features**
1. **Multi-Gene Networks**: >2 genes with complex dependencies
2. **Stochastic PRC2 Fluctuations**: Noise-driven boundary crossings
3. **Experimental Validation Framework**: Testable predictions for optogenetic control

## Expected Biological Insights

### **Boundary Crossing Scenarios**
1. **Low PRC2**: Dll4 (sustained) activates first, then Dll1 (pulsatile)
2. **High PRC2**: Dll1 (pulsatile) activates first, then Dll4 (sustained)  
3. **Critical PRC2**: Simultaneous activation or oscillatory switching

### **Sequential Activation Patterns**
1. **Feed-forward activation**: Gene 1 → enhances Gene 2 activation
2. **Feed-back inhibition**: Gene 2 → suppresses Gene 1 activation
3. **Bistable switching**: PRC2-dependent dominance patterns

This enhanced model would capture the **non-monotonic, boundary-crossing behavior** needed to understand how PRC2 levels control **sequential gene activation ordering** in developmental contexts.