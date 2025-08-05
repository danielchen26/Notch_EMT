using Catalyst, DataFrames, CSV

# Import the model
model_pulsatile = @reaction_network begin
    (A * (1 + sign(cos(w * t + ϕ))), 1.0), R ↔ NR               # NICD binds RBPJ
    (k1, k2), M + R ↔ MR          # MITF binds RBPJ
    k0, MR --> MR + KDM5A            # MITF-RBPJ recruit KDM5A
    d, H4 + KDM5A --> H0 + KDM5A       # Demethylation of active mark
    m, H0 + PRC2 --> H27 + PRC2   # Methylation to get repressive mark
    1.0, H27 + KDM6A --> H0 + KDM6A  # Demethylation of repressive mark
    1.0, H0 + KMT --> H4 + KMT    # Methylation to get active mark
    p, H27 --> H27 + PRC2           # PRC2 is enhenced by H27 mark
    kk, H4 --> H4 + KDM6A        # KDM6A is enhenced by H4 mark
    pp, H4 --> H4 + KMT          # KMT is enhenced by H4 mark
    k, H27 --> H27 + KDM5A                # KDM5A is enhenced by H27 mark
    δ, (PRC2, KDM5A, KDM6A, KMT) --> ∅                    # Degradation of histone reader and writers
    α1, ∅ --> (KDM6A, KMT, PRC2, KDM5A)
end

println("Species from model:")
for (i, s) in enumerate(species(model_pulsatile))
    println("  $i: $s")
end

println("\nParameters from model:")
for (i, p) in enumerate(parameters(model_pulsatile))
    println("  $i: $p")
end

# Load the database to check column names
db = CSV.File("../Notch_EMT_data/Notch_params_complete.csv") |> DataFrame
println("\nInitial condition columns from database:")
initi_names = names(db)[13:end]
for (i, name) in enumerate(initi_names)
    println("  $i: $name")
end

# Check db_idx 49 specifically
println("\nRow 49 from database:")
println("Parameters:")
for name in ["k0", "k1", "k2", "d", "m", "p", "k", "pp", "kk", "δ", "α1"]
    println("  $name = $(db[49, name])")
end

println("\nInitial conditions:")
for name in initi_names
    println("  $name = $(db[49, name])")
end

# Let's also check if prc2=0.41 makes sense
println("\nPRC2 parameter 'p' = $(db[49, :p])")
println("Note: In figure 3, prc2=0.41 is used, which likely overrides parameter 'p' or 'm' or both")