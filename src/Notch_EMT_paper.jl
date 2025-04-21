# src/Notch_EMT_paper.jl
# ===================================================================
# This script generates Figures 3, 4, and 5 from the Notch-EMT paper
# ===================================================================

## Load packages and setup environment
using Revise
using DifferentialEquations
using Plots; pyplot() # Use PyPlot consistently throughout
using Catalyst
using Catalyst: parameters
using DataFrames, CSV, Random, Distributions
using StatsPlots
using ProgressMeter
using Latexify, Measures, LaTeXStrings
includet("./Functions.jl")

## Initialize data and directories
println("Loading database and models...")
# Load database and parameters
db, p_names, initi_names, parameter_set, initial_condition, rand_idx_set, db_idx = loading_database()

# Create output directories
figure_path = joinpath(dirname(@__DIR__), "figures", "paper")
data_path = joinpath(dirname(@__DIR__), "Data", "paper_data")
isdir(figure_path) || mkpath(figure_path)
isdir(data_path) || mkpath(data_path)

# Import models
model_pulsatile = Import_model(; type="pulsatile")

## Figure 3: Sustained vs. Pulsatile Dynamics
function generate_figure3()
    println("Generating Figure 3: Sustained vs. Pulsatile dynamics...")
    
    # Common parameters
    db_idx = 49
    amplitude = 50.0
    T_init = 1e-10
    ΔT = 100
    prc2 = 0.41
    T_final = 1.5 * ΔT
    
    # Panel A: Sustained input (freq = 0.0)
    plt_Dll4 = single_solve_plot(; 
        model=model_pulsatile, 
        db_idx=db_idx, 
        phase=0, 
        freq=0.0, 
        prc2=prc2,
        amplitude=amplitude, 
        T_init=T_init, 
        ΔT=ΔT, 
        type="sustained", 
        T_final=T_final
    )
    
    # Improve panel A plot styling and select specific variables
    plot!(plt_Dll4,
        title="A=50.0, freq=0.0, ST=10.3",
        xlabel="Time",
        ylabel="Concentration",
        idxs=[1, 6, 9], # MR, H4, H27
        labels=["MR(t)" "H4(t)" "H27(t)"],
        legendfontsize=10,
        guidefontsize=12,
        linewidth=2
    )
    
    # Panel B: Pulsatile input (freq = 0.43)
    plt_Dll1 = single_solve_plot(; 
        model=model_pulsatile, 
        db_idx=db_idx, 
        phase=0, 
        freq=0.43, 
        prc2=prc2,
        amplitude=amplitude, 
        T_init=T_init, 
        ΔT=ΔT, 
        type="pulsatile", 
        T_final=T_final
    )
    
    # Improve panel B plot styling and select specific variables
    plot!(plt_Dll1,
        title="A=50.0, freq=0.43, ST=56.17",
        xlabel="Time",
        ylabel="Concentration",
        idxs=[1, 6, 9], # MR, H4, H27
        labels=["MR(t)" "H4(t)" "H27(t)"],
        legendfontsize=10,
        guidefontsize=12,
        linewidth=2
    )
    
    # Combine plots
    fig3 = plot(plt_Dll4, plt_Dll1, layout=(1,2), size=(1000, 400))
    
    # Save figure
    savefig(fig3, joinpath(figure_path, "figure3_sustained_vs_pulsatile.png"))
    println("Figure 3 saved.")
    
    return fig3
end

## Figure 4: A-ω Relation (Switching Boundary and ST-ω)
function generate_figure4()
    println("Generating Figure 4: A-ω relation...")
    
    # Parameters
    db_idx = 49  # Gene ID used in the paper
    amplitude_range = 0:1:300 # Use fine range for smooth plot
    freq_range = 0:0.02:2
    ΔT = 100
    
    # Define save paths using pathgen
    df_save_path, figure_save_path = pathgen(db_idx=db_idx, type="pulsatile")
    
    # Construct filename based on current ranges (like original)
    filename = joinpath(df_save_path, "freq :$freq_range" * "_|_" * "amplitude :$amplitude_range.csv")
    
    # Load or generate data
    local df
    if isfile(filename)
        println("Loading pre-computed data from: $(filename)")
        df = CSV.read(filename, DataFrame)
    else
        println("Generating data for Figure 4 (this may take time)...")
        df = A_ω_st_relation(;
            model=model_pulsatile,
            db_idx=db_idx,
            amplitude_range=amplitude_range,
            freq_range=freq_range,
            ΔT=ΔT
        )
        # Save data to the constructed path
        CSV.write(filename, df)
        println("Data saved to: $(filename)")
    end
    
    # Panel A: A-ω Switching Boundary
    println("Generating Panel 4A: A-ω Boundary")
    df_min_amp = extract_min_amp(df)
    
    fig4a = plot(df_min_amp.amp, df_min_amp.freq, 
        seriestype=:line,  
        linewidth=3,
        color=:red,
        legend=:topleft,
        label="Switching Boundary",
        xlabel=L"Driving Amplitude ($A$)", 
        ylabel=L"Driving Frequency ($\omega$)",
        guidefontsize=12,
        dpi=500
    )
    
    # Panel B: Frequency vs. Switching Time grouped by amplitude
    println("Generating Panel 4B: Frequency vs. Switching Time")
    fig4b = df_freq_vs_ST_groupby_amp(df; 
        amplitude_select=collect(50:50:300),
        figure_save_path=figure_save_path # Pass the path if needed, remove other style args
    )
    
    # Combine plots
    fig4 = plot(fig4a, fig4b, layout=(1,2), size=(1000, 400))
    
    # Save figure
    savefig(fig4, joinpath(figure_path, "figure4_A_omega_relation.png"))
    println("Figure 4 saved.")
    
    return fig4
end

## Figure 5: A-ω Curve Controlled by PRC2 Rate
function generate_figure5()
    println("Generating Figure 5: A-ω curves controlled by PRC2 rate...")
    
    # Parameters
    db_idx = 592 # Although data is loaded directly, keep params for potential generation
    amplitude_range = 0:5:500
    freq_range = 0:0.05:1
    prc2_range = 0.1:0.1:1
    ΔT = 100
    
    # Define the specific filename used in the original script
    target_filename = "df_592_3d.csv"
    # Look for the file in the main project directory (since it's attached there)
    data_file_path = joinpath(dirname(@__DIR__), target_filename)
    
    # Load or generate data
    local df_592_3d
    if isfile(data_file_path)
        println("Loading pre-computed data from: $(data_file_path)")
        df_592_3d = CSV.read(data_file_path, DataFrame)
    else
        println("Data file $(data_file_path) not found.")
        println("Generating data for Figure 5 (this may take significant time)...")
        df_592_3d = A_ω_st_relation_prc2_range(;
            model=model_pulsatile,
            db_idx=db_idx,
            amplitude_range=amplitude_range,
            freq_range=freq_range,
            prc2_range=prc2_range,
            ΔT=ΔT
        )
        # Save the newly generated data to the main Data path for future use
        CSV.write(data_file_path, df_592_3d)
        println("Generated data saved to: $(data_file_path)")
    end
    
    fig5 = plot_all_prc2(df_592_3d, :prc2, :amp, :freq; 
        legend=:outerright, 
        foreground_color_legend=nothing, 
        legendfontsize=9,
        legendtitle="PRC2 Rate",
        guidefontsize=12,
        xlabel="Driving Amplitude (A)",
        ylabel=L"Driving Frequency ($\omega$)",
        title="A-ω Decision Boundary Regulated by PRC2 Rate",
        dpi=500,
        size=(650, 450)
    )
    
    # Save figure
    savefig(fig5, joinpath(figure_path, "figure5_A_omega_vs_PRC2.png"))
    println("Figure 5 saved.")
    
    return fig5
end


## Generate just Figure 3
generate_figure3()

## Generate just Figure 4
generate_figure4()

## Generate just Figure 5
generate_figure5() 
