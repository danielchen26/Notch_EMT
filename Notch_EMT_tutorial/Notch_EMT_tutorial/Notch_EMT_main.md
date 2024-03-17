```@meta
EditURL = "../Notch_EMT_main.jl"
```

````@example Notch_EMT_main
# ======== Loading Pkgs and functions =========
@time include("functions/preload.jl")
show_avail_functions()
# =============================================

# ======= Loading model ======
@time include("functions/models.jl")
# =================================
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

