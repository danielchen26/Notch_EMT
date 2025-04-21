# src/Notch_EMT_paper.jl
# ===================================================================
# This script generates Figures 3, 4, and 5 from the Notch-EMT paper
# by simulating Notch signaling dynamics under different conditions.
# It utilizes pre-computed data when available to speed up figure generation.
# ===================================================================

# -------------------------------------------------------------------
# SECTION 1: Load Packages and Setup Environment
# -------------------------------------------------------------------
# Load necessary Julia packages for differential equations, plotting,
# data manipulation, and system modeling.
using Revise # Helps with interactive code development
using DifferentialEquations # For solving ODEs
using Plots; pyplot() # Plotting library; set PyPlot as the backend for consistency
using Catalyst # For defining reaction networks
using Catalyst: parameters # Accessing model parameters
using DataFrames, CSV, Random, Distributions # Data handling and statistics
using StatsPlots # Statistical plotting recipes
using ProgressMeter # Displaying progress bars for long computations
using Latexify, Measures, LaTeXStrings # For LaTeX rendering in plots

# Include custom functions defined in Functions.jl (e.g., model loading, simulation helpers, plotting)
includet("./Functions.jl")

# -------------------------------------------------------------------
# SECTION 2: Initialize Data, Directories, and Models
# -------------------------------------------------------------------
println("Loading database and models...")
# Load the parameter database, initial conditions, and get a random index
# (Note: the random db_idx is global, figure functions use local specific indices)
db, p_names, initi_names, parameter_set, initial_condition, rand_idx_set, global_db_idx = loading_database()

# Define paths for saving figures and potentially loading/saving data
figure_path = joinpath(dirname(@__DIR__), "figures", "paper") # Path for final paper figures
data_path_regular = joinpath(dirname(@__DIR__), "Data", "regular") # Path for general simulation data
data_path_main = joinpath(dirname(@__DIR__), "Data") # Main data directory

# Create output directories if they don't exist
isdir(figure_path) || mkpath(figure_path)
isdir(data_path_regular) || mkpath(data_path_regular) # Ensure regular data path exists
isdir(data_path_main) || mkpath(data_path_main) # Ensure main data path exists

# Import the Catalyst model for pulsatile Notch signaling
model_pulsatile = Import_model(; type="pulsatile")

# -------------------------------------------------------------------
# SECTION 3: Figure Generation Functions
# -------------------------------------------------------------------

# =======================
# Figure 3: Sustained vs. Pulsatile Dynamics
# =======================
# Compares the system's response to a sustained vs. a specific pulsatile input.
function generate_figure3()
    println("Generating Figure 3: Sustained vs. Pulsatile dynamics...")
    
    # Define common parameters for both panels
    db_idx = 49      # Specific gene parameter set ID
    amplitude = 50.0 # Input signal amplitude
    T_init = 1e-10   # Time before signal starts (must be > 0)
    ΔT = 100         # Duration of the signal pulse
    prc2 = 0.41      # Specific PRC2 rate for this simulation
    T_final = 1.5 * ΔT # Total simulation time relative to pulse duration
    
    # Ensure PyPlot backend is active for consistency
    pyplot()
    
    # --- Panel A: Sustained Input --- 
    # Simulate with frequency = 0.0 (sustained signal)
    plt_Dll4 = single_solve_plot(; 
        model=model_pulsatile, 
        db_idx=db_idx, 
        phase=0, 
        freq=0.0,        # Key parameter for sustained input
        prc2=prc2,
        amplitude=amplitude, 
        T_init=T_init, 
        ΔT=ΔT, 
        type="sustained", # Specifies signal type for plotting
        T_final=T_final
    )
    
    # Customize Panel A plot appearance
    plot!(plt_Dll4,
        title="A=50.0, freq=0.0, ST=10.3", # Title includes key parameters/results
        xlabel="Time",
        ylabel="Concentration",
        idxs=[1, 6, 9], # Select specific species indices to plot (MR, H4, H27)
        labels=["MR(t)" "H4(t)" "H27(t)"], # Labels for selected species
        legendfontsize=10,
        guidefontsize=14,
        tickfontsize=12,
        linewidth=2
    )
    # Save Panel A individually
    savefig(plt_Dll4, joinpath(figure_path, "figure3_panelA_sustained.png"))
    println("Figure 3 Panel A saved.")

    # --- Panel B: Pulsatile Input --- 
    # Simulate with a specific frequency (0.43)
    plt_Dll1 = single_solve_plot(; 
        model=model_pulsatile, 
        db_idx=db_idx, 
        phase=0, 
        freq=0.43,       # Key parameter for pulsatile input
        prc2=prc2,
        amplitude=amplitude, 
        T_init=T_init, 
        ΔT=ΔT, 
        type="pulsatile",  # Specifies signal type for plotting
        T_final=T_final
    )
    
    # Customize Panel B plot appearance
    plot!(plt_Dll1,
        title="A=50.0, freq=0.43, ST=56.17", # Title includes key parameters/results
        xlabel="Time",
        ylabel="Concentration",
        idxs=[1, 6, 9], # Select specific species indices to plot (MR, H4, H27)
        labels=["MR(t)" "H4(t)" "H27(t)"], # Labels for selected species
        legendfontsize=10,
        guidefontsize=14,
        tickfontsize=12,
        linewidth=2
    )
    # Save Panel B individually
    savefig(plt_Dll1, joinpath(figure_path, "figure3_panelB_pulsatile.png"))
    println("Figure 3 Panel B saved.")
    
    # --- Combine and Save --- 
    # Arrange Panel A and B side-by-side
    fig3 = plot(plt_Dll4, plt_Dll1, layout=(1,2), size=(1000, 400))
    
    # Save the combined figure to the designated path
    savefig(fig3, joinpath(figure_path, "figure3_sustained_vs_pulsatile.png"))
    println("Figure 3 combined saved.")
    
    return fig3 # Return the plot object
end

# =======================
# Figure 4: A-ω Relation (Switching Boundary and ST-ω)
# =======================
# Explores how switching behavior depends on signal amplitude (A) and frequency (ω).
function generate_figure4()
    println("Generating Figure 4: A-ω relation...")
    
    # Define parameters for the simulation sweep
    db_idx = 49      # Specific gene parameter set ID
    amplitude_range = 0:1:300 # Fine amplitude range for smooth boundary plot
    freq_range = 0:0.02:2   # Frequency range to test
    ΔT = 100         # Signal pulse duration
    
    # Define paths using the standard pathgen function (saves to Data/regular/)
    df_save_path, figure_save_path = pathgen(db_idx=db_idx, type="pulsatile")
    
    # Construct the expected filename based on current parameter ranges
    filename = joinpath(df_save_path, "freq :$freq_range" * "_|_" * "amplitude :$amplitude_range.csv")
    
    # Load data if it exists, otherwise generate it
    local df # Dataframe to hold A-ω-ST results
    if isfile(filename)
        println("Loading pre-computed data from: $(filename)")
        df = CSV.read(filename, DataFrame)
    else
        println("Generating data for Figure 4 (this may take time due to fine amplitude range)...")
        # Perform simulations across the amplitude and frequency ranges
        df = A_ω_st_relation(;
            model=model_pulsatile,
            db_idx=db_idx,
            amplitude_range=amplitude_range,
            freq_range=freq_range,
            ΔT=ΔT
        )
        # Save the newly generated data
        CSV.write(filename, df)
        println("Data saved to: $(filename)")
    end
    
    # Ensure PyPlot backend is active
    pyplot()

    # --- Panel A: A-ω Switching Boundary --- 
    println("Generating Panel 4A: A-ω Boundary")
    df_min_amp = extract_min_amp(df)
    
    # Create the base plot object
    fig4a = plot(
        xlabel=L"Driving Amplitude ($A$)", 
        ylabel=L"Driving Frequency ($\omega$)",
        guidefontsize=14, 
        tickfontsize=12,  
        legend=:topleft,
        dpi=500
    )

    # Plot the boundary line
    plot!(fig4a, df_min_amp.amp, df_min_amp.freq, 
        seriestype=:line,  
        linewidth=3,
        color=:red,
        label="Switching Boundary"
    )
    
    # Get plot limits for filling and annotation placement
    xlims_vals = xlims(fig4a)
    ylims_vals = ylims(fig4a) 

    # Fill above the line - visually distinct region
    plot!(fig4a, df_min_amp.amp, df_min_amp.freq, 
        fillrange = ylims_vals[2], 
        fillalpha = 0.25, 
        fillcolor = :lightblue, # Color for visual separation
        seriestype = :line, 
        linewidth = 0, 
        label = "" 
    )

    # Fill below the line - visually distinct region
    plot!(fig4a, df_min_amp.amp, df_min_amp.freq, 
        fillrange = 0, 
        fillalpha = 0.25, 
        fillcolor = :lightsalmon, # Color for visual separation
        seriestype = :line, 
        linewidth = 0, 
        label = "" 
    )

    # Add annotations based on correct logical regions (Left/Right)
    # Place "Non-Switching" clearly to the LEFT of the curve and in the upper (blue) region
    annotate!(fig4a, 
        xlims_vals[1] + 0.20 * (xlims_vals[2] - xlims_vals[1]), # x-position (left side, 20% - slightly adjusted)
        ylims_vals[1] + 0.80 * (ylims_vals[2] - ylims_vals[1]),  # y-position (upper part, 80%)
        text("Non-Switching Region", 12, :darkorange, :center) # Reverted to standard, size 12
    )
    # Place "Switching" clearly to the RIGHT of the curve and in the lower-left (orange) region
    annotate!(fig4a, 
        xlims_vals[1] + 0.65 * (xlims_vals[2] - xlims_vals[1]), # x-position (moved left, 65%)
        ylims_vals[1] + 0.20 * (ylims_vals[2] - ylims_vals[1]),  # y-position (lower part, 20%)
        text("Switching Region", 12, :darkblue, :center) # Reverted to standard, size 12
    )

    # Save Panel A individually
    savefig(fig4a, joinpath(figure_path, "figure4_panelA_boundary.png"))
    println("Figure 4 Panel A saved.")
    
    # --- Panel B: Frequency vs. Switching Time (ST) --- 
    println("Generating Panel 4B: Frequency vs. Switching Time")
    # Plot ST vs. frequency, grouped by selected amplitudes
    fig4b = df_freq_vs_ST_groupby_amp(df; 
        amplitude_select=collect(50:50:300), # Select specific amplitudes to show
        figure_save_path=figure_save_path # Pass path for potential internal saving
    )
    # Save Panel B individually
    # Note: df_freq_vs_ST_groupby_amp might already save this if figure_save_path is provided
    # Adding an explicit save here ensures it's saved with a specific name.
    savefig(fig4b, joinpath(figure_path, "figure4_panelB_ST_vs_omega.png"))
    println("Figure 4 Panel B saved.")
    
    # --- Combine and Save --- 
    # Arrange Panel A and B side-by-side
    fig4 = plot(fig4a, fig4b, layout=(1,2), size=(1000, 400))
    
    # Save the combined figure
    savefig(fig4, joinpath(figure_path, "figure4_A_omega_relation.png"))
    println("Figure 4 combined saved.")
    
    return fig4 # Return the plot object
end

# =======================
# Figure 5: A-ω Curve Controlled by PRC2 Rate
# =======================
# Investigates how the A-ω switching boundary changes with varying PRC2 rates.
function generate_figure5()
    println("Generating Figure 5: A-ω curves controlled by PRC2 rate...")
    
    # Parameters
    db_idx = 592 # Although data is loaded directly, keep params for potential generation
    amplitude_range = 0:5:500
    freq_range = 0:0.05:1
    prc2_range = 0.1:0.1:1 # Range of PRC2 rates to test
    ΔT = 100
    
    # Define the specific filename used in the original script
    target_filename = "df_592_3d.csv"
    # Ensure we look for the file in the main project directory (project root)
    data_file_path = joinpath(dirname(@__DIR__), target_filename) # Correct path
    
    # Load the specific curated data file if it exists, otherwise generate it.
    local df_592_3d # Dataframe holding A-ω-ST results across PRC2 rates
    if isfile(data_file_path)
        println("Loading pre-computed data from: $(data_file_path)")
        df_592_3d = CSV.read(data_file_path, DataFrame)
    else
        println("Data file $(data_file_path) not found.") # This message uses the correct path
        println("Generating data for Figure 5 (this may take significant time)...")
        # Perform simulations across amplitude, frequency, and PRC2 ranges
        df_592_3d = A_ω_st_relation_prc2_range(;
            model=model_pulsatile,
            db_idx=db_idx,
            amplitude_range=amplitude_range,
            freq_range=freq_range,
            prc2_range=prc2_range,
            ΔT=ΔT
        )
        # Save the newly generated data to the main project path for future use
        CSV.write(data_file_path, df_592_3d)
        println("Generated data saved to: $(data_file_path)")
    end
    
    # Ensure PyPlot backend is active
    pyplot()
    
    # --- Generate the Plot --- 
    # Use the dedicated plotting function (defined in Functions.jl)
    fig5 = plot_all_prc2(df_592_3d, :prc2, :amp, :freq; # Plot amp vs freq, grouped by prc2 
        legend=:outerright, 
        foreground_color_legend=nothing, 
        legendfontsize=9,
        legendtitle="PRC2 Rate",
        guidefontsize=14,
        tickfontsize=12,
        xlabel="Driving Amplitude (A)",
        ylabel=L"Driving Frequency ($\omega$)",
        title="A-ω Decision Boundary Regulated by PRC2 Rate",
        dpi=500,
        size=(650, 450)
    )
    
    # --- Save Figure --- 
    savefig(fig5, joinpath(figure_path, "figure5_A_omega_vs_PRC2.png"))
    println("Figure 5 saved.")
    
    return fig5 # Return the plot object
end

# -------------------------------------------------------------------
# SECTION 4: Script Execution Logic
# -------------------------------------------------------------------

# Define a function to run all figure generations sequentially (useful for batch processing)
function run_all()
    println("Generating all figures...")
    fig3 = generate_figure3()
    fig4 = generate_figure4()
    fig5 = generate_figure5()
    println("All figures generated and saved.")
    return fig3, fig4, fig5 # Return all plot objects
end

# --- How to Use This Script --- 
# Option 1: Interactive Use (Julia REPL)
#   - Start Julia: `julia`
#   - Include this file: `include("src/Notch_EMT_paper.jl")`
#   - Call functions individually: `generate_figure3()`, `generate_figure4()`, `generate_figure5()`, or `run_all()`
#
# Option 2: Direct Execution (Command Line or Editor)
#   - Uncomment ONE of the lines below to generate a specific figure when running the script directly.
#   - Example: `julia src/Notch_EMT_paper.jl`

# --- Uncomment ONE line below to run directly --- 
# generate_figure3()
# generate_figure4()
# generate_figure5() 
run_all() # Uncomment this to generate all figures

# Print guidance message when the script is included or run without uncommenting a generation line
# Simplified condition to check if running as main script non-interactively
if !isinteractive() && abspath(PROGRAM_FILE) == @__FILE__
    println("\nScript finished setup. Check comments at the end for how to generate figures.")
    println("To generate figures:")
    println("  1. Uncomment a line like 'generate_figure3()' or 'run_all()' at the end of this script and run again.")
    println("  2. Or, run interactively: include(\"src/Notch_EMT_paper.jl\") then call functions like generate_figure3().")
end 
