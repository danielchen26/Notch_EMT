```@meta
EditURL = "Notch_EMT_main.jl"
```

````@example Notch_EMT_main
# ======== Loading Pkgs and functions =========
@time include(joinpath(@__DIR__, "functions/preload.jl"))
include(joinpath(@__DIR__, "functions/helper_fun.jl"))
show_avail_functions()
# =============================================




# ======= Loading model ======
@time include(joinpath(@__DIR__, "functions/models.jl"))
# =================================




# ======= Single gene case =======
run_path = joinpath(@__DIR__, "functions/single_gene_case.jl")
include(run_path);
run_and_check_new_vars(run_path)
# ============================================================
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

