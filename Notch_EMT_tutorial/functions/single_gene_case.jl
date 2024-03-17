using Logging: @info

function process_signal(db_idx, freq, phase, amplitude, T_init, ΔT, tspan, model, db, p_names)
    @info "Processing signal" db_idx freq phase amplitude T_init ΔT tspan
    @info "Model selected: $(latexify(model))"
    u0map, pmap, p, tspan, ts, cb, sol = single_solve(; db_idx=db_idx, freq=freq, phase=phase, amplitude=amplitude, T_init=T_init, ΔT=ΔT, tspan=tspan, phase_reset=true);
    @info "Switching time:" t_switching = find_ST(sol)
    plt = plot(sol, idxs=[4, 5, 6, 9], lw=1.5, xlabel="Time", ylabel="Concentration", dpi=500)
    osci_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)));
    tt = ts[1]:0.01:ts[2];
    osci_signal.(tt, amplitude, freq, 0.0)
    plot!(plt, tt, osci_signal.(tt, amplitude, freq, 0.0),
        label="Pulsatile Input", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
        fill=(0, 0.3, :darkgreen), color="black", dpi=300)
    display(plt)
    return plt, sol
end


## ======== single gene case for id : 592 =====
db_idx = 592
freq = 0.0
phase = 0
amplitude = 220
T_init = 1e-10
ΔT = 100
tspan = (0.0, 350.0)
model = model_pulsatile

# An example instance after the function
example_instance = process_signal(db_idx, freq, phase, amplitude, T_init, ΔT, tspan, model, db, p_names); @nl


##💚 This is single gene model comparison (we are using pulsatile signal)
## ================== comparison of pulsatile and bump model =================
@info "Comparing pulsatile and bump model for single gene 49"
function S_gene_model_comparison(;db_idx=49, phase=0, freq=0.5, amplitude=65.0, T_init=1e-10, ΔT=80)
    # Pulsatile model
    @info "Pulsatile model"
    model = model_pulsatile
    signal_type = "pulsatile"
    plt_pulsatile = single_solve_plot(; model=model, db_idx=db_idx, phase=phase, freq=freq, amplitude=amplitude, T_init=T_init, ΔT=ΔT, type=signal_type)

    # Bump model
    @info "Bump model"
    model = model_bump
    signal_type = "bump"
    plt_bump = single_solve_plot(; model=model, db_idx=db_idx, phase=phase, freq=freq, amplitude=amplitude, T_init=T_init, ΔT=ΔT, type=signal_type)

    # Combine plots
    plt_2model = plot(plt_pulsatile, plt_bump, layout=(2, 1))
    display(plt_2model)
    return plt_2model
end

example_instance = S_gene_model_comparison(; db_idx=49, phase=0, freq=0.5, amplitude=65.0, T_init=1e-10, ΔT=80)
@nl

## ================== test discontinuity =================
function test_discontinuity(;ΔT_range = 0.01:10:1000, T_init = 0.01)
    #! Have to set T_init to avoid the discontinuity
    for ΔT ∈ ΔT_range
        plt = S_gene_model_comparison(ΔT=ΔT, T_init=T_init)
        display(plt)
    end
end

@info " `Test discontinuity` function available to see if the switching time is correct as ΔT increases and T_init changes" 
# test_discontinuity();


## 💚 This is single gene signal comparison, assuiming pusatile model.
## =========== Single gene , pulsatile model DLL1 vs DLL4 comparison ==========

"""
`S_gene_Dll1vsDll4` is a function to plot a single gene with ID 
default to 592 comparing sustained signal vs pulsatile signal.

# Parameters
- `db_idx`: The gene ID. Default is 592.
- `phase`: The phase of the signal. Default is 0.
- `freq`: The frequency of the signal. Default is 0.0.
- `amplitude`: The amplitude of the signal. Default is 165.0.
- `T_init`: The initial time. Default is 0.01.
- `ΔT`: The time interval. Default is 100.

# Returns
- A plot that shows the comparison between sustained and pulsatile signals of a single gene.
"""
function S_gene_Dll1vsDll4(; db_idx=592, phase=0, freq=0.0, amplitude=165.0, T_init=0.01, ΔT=100)
    plt_sustained = single_solve_plot(; 
        model=model_pulsatile, 
        db_idx=db_idx, 
        phase=0.0, 
        freq=0.0, 
        amplitude=amplitude, 
        T_init=T_init, 
        ΔT=ΔT, 
        type="sustained", 
        phase_reset=true
    )
    plt_pulsatile = single_solve_plot(; 
        model=model_pulsatile, 
        db_idx=db_idx, 
        phase=phase, 
        freq=freq, 
        amplitude=amplitude, 
        T_init=T_init, 
        ΔT=ΔT, 
        type="pulsatile", 
        phase_reset=true
    )
    return plot(plt_sustained, plt_pulsatile, layout=(2, 1))
end



#! need to add args.
function test_S_gene_Dll1vsDll4()
    # * showed some freq and amplitude allows pulsatile signal switch states after signal is off (between pulses)
    for freq ∈ 0.01:0.05:1#, amplitude ∈ 165:10:300
        plt = S_gene_Dll1vsDll4(ΔT=100, T_init=0.01, freq=freq, amplitude=300)#, amplitude = amplitude)
        display(plt)
        sleep(0.1)
    end
end
# using Test
# @test test_S_gene_Dll1vsDll4()



## =========== Single gene, pulsatile model DLL1 vs DLL4 comparison ==========
## ==== Dll4 and Dll1 can both induce sustained histone state ====
#🔴 Figure 5 a, b
function figure5(model, T_init, signal_type; savefig=true)
    plt_Dll1 = single_solve_plot(; 
        model=model, 
        db_idx=49, 
        phase=0, 
        freq=0.43, 
        prc2=0.41, 
        amplitude=50.0, 
        T_init=T_init, 
        ΔT=100, 
        type=signal_type, 
        T_final=1.5 * ΔT
    )

    plt_Dll4 = single_solve_plot(; 
        model=model, 
        db_idx=49, 
        phase=0, 
        freq=0.0, 
        prc2=0.41, 
        amplitude=50.0, 
        T_init=T_init, 
        ΔT=100, 
        type=signal_type, 
        T_final=1.5 * ΔT
    )

    if savefig
        savefig(plt_Dll1, "./figures/Dll1_multiple.png")
        savefig(plt_Dll4, "./figures/Dll4.png")
    else
        plot(plt_Dll4, plt_Dll1, layout=(1,2), size=(800, 400))
    end
end

figure5(model_pulsatile, 1e-10, "pulsatile", savefig=false)

