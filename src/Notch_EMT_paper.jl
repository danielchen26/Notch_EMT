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
db, p_names, initi_names, parameter_set, initial_condition, rand_idx_set, global_db_idx = loading_database(; data_path="../Notch_EMT_data/Notch_params_complete.csv")

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
    # Use the single_solve_plot wrapper function
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
        T_final=T_final,
        title="off" # Turn off default title in wrapper
    )
    
    # Apply styling and attempt to select specific variables using plot!
    # WARNING: This might not remove KDM5A from the legend/plot reliably.
    plot!(plt_Dll4,
        title=L"A=50.0, \omega=0.0, ST=10.3",
        xlabel=L"Time ($T$)",
        ylabel="Concentration",
        idxs=[4, 6, 9], # Attempt to select MR, H4, H27
        labels=["MR(t)" "H4(t)" "H27(t)"],
        legendfontsize=10,
        guidefontsize=18,  
        tickfontsize=14,   
        titlefontsize=18,  
        legend_framestyle=:none, 
        linewidth=2
    )
    
    # Save Panel A individually
    savefig(plt_Dll4, joinpath(figure_path, "figure3_panelA_sustained.png"))
    println("Figure 3 Panel A saved.")

    # --- Panel B: Pulsatile Input --- 
    # Use the single_solve_plot wrapper function
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
        T_final=T_final,
        phase_reset=true,
        title="off" # Turn off default title in wrapper
    )
    
    # Apply styling and attempt to select specific variables using plot!
    # WARNING: This might not remove KDM5A from the legend/plot reliably.
    plot!(plt_Dll1,
        title=L"A=50.0, \omega=0.43, ST=56.17",
        xlabel=L"Time ($T$)",
        ylabel="Concentration",
        idxs=[4, 6, 9], # Attempt to select MR, H4, H27
        labels=["MR(t)" "H4(t)" "H27(t)"],
        legendfontsize=10,
        guidefontsize=18,  
        tickfontsize=14,   
        titlefontsize=18,  
        legend_framestyle=:none, 
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
    
    # Find the minimum amplitude required for switching at frequency zero (sustained input)
    sustained_switching_amp_row = filter(row -> row.freq == 0, df_min_amp)
    sustained_amp = isempty(sustained_switching_amp_row) ? missing : sustained_switching_amp_row[1, :amp]
    max_amp_plot = maximum(df_min_amp.amp) # Get max amplitude for xlim padding

    # Create the base plot object
    fig4a = plot(
        xlabel=L"Driving Amplitude ($A$)", 
        ylabel=L"Driving Frequency ($\omega$)",
        guidefontsize=18,
        tickfontsize=14,
        legendfontsize=10,
        legend=:topleft,
        xlims=(0, max_amp_plot * 1.05), # Ensure x-axis starts at 0 and add 5% padding
        dpi=500
    )
    
    # Plot the boundary line
    plot!(fig4a, df_min_amp.amp, df_min_amp.freq, 
        seriestype=:line,  
        linewidth=3,
        color=:red,
        label="Switching Boundary"
    )
    
    # Add point and annotation for sustained switching amplitude if found
    if !ismissing(sustained_amp)
        scatter!(fig4a, [sustained_amp], [0], 
            markercolor=:green,
            markersize=5, 
            label="Min. Sustained Amp")
        # Slightly offset annotation for clarity
        annotate!(fig4a, sustained_amp + 0.02 * xlims(fig4a)[2], 0.05 * ylims(fig4a)[2], # Increased x-offset slightly
            text("($sustained_amp, 0)", 10, :green, :left)) # Changed color to :green
    end

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
    # Look for the file in the Data/precomputed directory
    data_file_path = joinpath(dirname(@__DIR__), "Data", "precomputed", target_filename)
    
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
        legendfontsize=10, # Reduced from 12
        legendtitle="PRC2 Rate",
        guidefontsize=18,
        tickfontsize=14,
        xlabel=L"Driving Amplitude ($A$)",
        ylabel=L"Driving Frequency ($\omega$)",
        title=L"$A$-$\omega$ Decision Boundary Regulated by PRC2 Rate",
        titlefontsize=18,
        dpi=500,
        size=(650, 450)
    )
    
    # --- Save Figure --- 
    savefig(fig5, joinpath(figure_path, "figure5_A_omega_vs_PRC2.png"))
    println("Figure 5 saved.")
    
    return fig5 # Return the plot object
end

# -------------------------------------------------------------------
# SECTION 3.5: Custom Simulation for Critical Amplitude Check
# -------------------------------------------------------------------

# =======================
# Check Critical Amplitude for Specific Parameters
# =======================
# Performs a simulation to find the critical amplitude for a given
# db_idx, target frequency, and new ΔT.
function check_critical_amplitude_custom(; model=model_pulsatile, db_idx=49, target_freq=1.0, new_ΔT=200.0, amplitude_range=0:1:300, prc2_val="NA")
    println("\nChecking critical amplitude for: db_idx=$db_idx, ω=$target_freq, ΔT=$new_ΔT, PRC2=$prc2_val")

    # Generate A-ω-ST data for the specific frequency and ΔT
    # A_ω_st_relation returns a DataFrame of *switching* events
    df_switching_events = A_ω_st_relation(;
        model=model,
        db_idx=db_idx,
        amplitude_range=amplitude_range,
        freq_range=[target_freq], # Pass the target frequency as a single-element array
        ΔT=new_ΔT,
        prc2=prc2_val, # Use the prc2 value from Figure 3 & 4 if not specified, or allow override
        mute_parameter_disp=true
    )

    if isempty(df_switching_events)
        println("No switching events found for the given parameters and amplitude range.")
        println("Critical amplitude might be higher than $(maximum(amplitude_range)) or switching does not occur under these conditions.")
        return nothing
    end

    # Filter for the target frequency (should be redundant if freq_range was [target_freq] but good for safety)
    df_target_freq = filter(row -> row.freq == target_freq, df_switching_events)

    if isempty(df_target_freq)
        println("No switching events found specifically for frequency ω=$target_freq after filtering.")
        # This case should ideally not be reached if A_ω_st_relation works as expected with a single freq in freq_range
        return nothing
    end

    # The critical amplitude is the minimum amplitude in this filtered DataFrame
    critical_amplitude = minimum(df_target_freq.amp)

    println("Critical amplitude for ω=$target_freq and ΔT=$new_ΔT (db_idx=$db_idx, PRC2=$prc2_val): $critical_amplitude")
    
    # For comparison, let's find the critical amplitude for the original ΔT=100 at the same frequency
    println("\nFor comparison, checking critical amplitude for original ΔT=100 at ω=$target_freq (db_idx=$db_idx, PRC2=$prc2_val)")
    df_switching_events_original_ΔT = A_ω_st_relation(;
        model=model,
        db_idx=db_idx,
        amplitude_range=amplitude_range,
        freq_range=[target_freq],
        ΔT=100.0, # Original ΔT from Figure 4
        prc2=prc2_val,
        mute_parameter_disp=true
    )

    if isempty(df_switching_events_original_ΔT)
        println("No switching events found for original ΔT=100 at ω=$target_freq.")
    else
        df_target_freq_original_ΔT = filter(row -> row.freq == target_freq, df_switching_events_original_ΔT)
        if isempty(df_target_freq_original_ΔT)
            println("No switching events found specifically for frequency ω=$target_freq with original ΔT=100 after filtering.")
        else
            critical_amplitude_original_ΔT = minimum(df_target_freq_original_ΔT.amp)
            println("Critical amplitude for ω=$target_freq and original ΔT=100 (db_idx=$db_idx, PRC2=$prc2_val): $critical_amplitude_original_ΔT")
        end
    end
    
    return critical_amplitude
end

# =======================
# Helper Function to Get A-ω Boundary Data for a specific ΔT
# =======================
# Generates or loads A-ω boundary data for a given db_idx, amplitude range, 
# frequency range, ΔT, and prc2 value.
function get_A_omega_boundary_data(; model=model_pulsatile, db_idx=49, amplitude_range=0:1:300, freq_range=0:0.02:2, ΔT_value=100.0, prc2_val="NA")
    println("\nGenerating or loading A-ω boundary data for: db_idx=$db_idx, ΔT=$ΔT_value, PRC2=$prc2_val")
    println("Amplitude range: $amplitude_range, Frequency range: $freq_range")

    # Define paths using the standard pathgen function (saves to Data/regular/db_idx:X/)
    df_save_path, _ = pathgen(db_idx=db_idx, type="pulsatile_deltaT_effect") # Using a slightly different type to avoid filename clashes if needed

    # Sanitize range string representations for filename (replace : with -)
    freq_range_str =replace(string(freq_range), ":" => "-")
    amp_range_str = replace(string(amplitude_range), ":" => "-")
    
    # Construct the expected filename based on parameters
    filename = joinpath(df_save_path, "A_omega_boundary_db_idx_$(db_idx)_freq_$(freq_range_str)_amp_$(amp_range_str)_deltaT_$(ΔT_value)_prc2_$(prc2_val).csv")
    
    local df_boundary # Dataframe to hold A-ω boundary
    if isfile(filename)
        println("Loading pre-computed A-ω boundary data from: $(filename)")
        df_boundary = CSV.read(filename, DataFrame)
    else
        println("Generating A-ω switching data (this may take time)...")
        # Perform simulations across the amplitude and frequency ranges for the given ΔT
        df_switching_events = A_ω_st_relation(
            model=model,
            db_idx=db_idx,
            amplitude_range=amplitude_range,
            freq_range=freq_range,
            ΔT=ΔT_value,
            prc2=prc2_val,
            mute_parameter_disp=true
        )
        
        if isempty(df_switching_events)
            println("Warning: No switching events found for db_idx=$db_idx, ΔT=$ΔT_value, PRC2=$prc2_val with the given ranges.")
            println("Cannot generate A-ω boundary. Returning an empty DataFrame.")
            # Return an empty DataFrame with expected columns if no switching occurs
            return DataFrame(amp=Float64[], freq=Float64[])
        end

        println("Extracting minimum amplitude for switching boundary...")
        df_boundary = extract_min_amp(df_switching_events)
        
        # Save the newly generated boundary data
        mkpath(dirname(filename)) # Ensure directory exists
        CSV.write(filename, df_boundary)
        println("A-ω boundary data saved to: $(filename)")
    end
    return df_boundary
end

# =======================
# Figure Compare A-ω Curves for Different ΔT
# =======================
# Generates and plots A-ω switching boundary curves for different ΔT values
# to visualize the effect of pulse duration on the switching threshold.
function generate_figure_compare_A_omega_for_different_deltaT()
    println("\nGenerating Figure: Comparison of A-ω boundaries for different ΔT values...")
    
    # Common parameters for the comparison
    db_idx = 49      # Specific gene parameter set ID (same as Fig 4)
    amplitude_range = 0:1:300 # FINER amplitude range (step 1, up to 300)
    freq_range = 0:0.02:2   # FINER frequency range (step 0.02)
    prc2_val = "NA"       # Use default PRC2 for db_idx 49, consistent with Fig 4
    
    # ΔT values to compare
    deltaT_values = [100.0, 200.0] # Original ΔT and the new one

    # Ensure PyPlot backend is active
    pyplot()
    
    # Create the base plot object
    fig_compare_deltaT = plot(
        xlabel=L"Driving Amplitude ($A$)", 
        ylabel=L"Driving Frequency ($\omega$)",
        title="", # REMOVED title
        guidefontsize=18,
        tickfontsize=14,
        legendfontsize=10,
        legend=:outerright, # Legend outside and centered on the right
        dpi=500
    )
    
    max_amp_overall = 0.0 # To adjust xlims later

    # Colors for different ΔT curves
    colors = [:blue, :red, :green, :purple] # Add more if comparing more ΔT values

    for (idx, deltaT) in enumerate(deltaT_values)
        println("\nProcessing for ΔT = $deltaT...")
        df_boundary = get_A_omega_boundary_data(
            db_idx=db_idx,
            amplitude_range=amplitude_range,
            freq_range=freq_range,
            ΔT_value=deltaT,
            prc2_val=prc2_val
        )
        
        if !isempty(df_boundary) && nrow(df_boundary) > 0
            # Sort by frequency for correct line plotting, though extract_min_amp usually does this
            sort!(df_boundary, :freq)
            
            plot!(fig_compare_deltaT, df_boundary.amp, df_boundary.freq, 
                seriestype=:line,  
                linewidth=2.5,
                color=colors[idx % length(colors) + 1], # Cycle through colors
                label="ΔT = $deltaT"
            )
            current_max_amp = maximum(df_boundary.amp)
            if current_max_amp > max_amp_overall
                max_amp_overall = current_max_amp
            end
        else
            println("No boundary data to plot for ΔT = $deltaT.")
        end
    end
    
    # Adjust xlims if data was plotted
    if max_amp_overall > 0
        xlims!(fig_compare_deltaT, (0, max_amp_overall * 1.05))
    else # Default xlims if no data
        xlims!(fig_compare_deltaT, (0, maximum(amplitude_range) * 1.05))
    end

    # Save the figure
    fig_save_path = joinpath(figure_path, "figure_compare_A_omega_deltaT_db$(db_idx).png")
    savefig(fig_compare_deltaT, fig_save_path)
    println("Comparison figure saved to: $fig_save_path")
    
    return fig_compare_deltaT
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

# Perform the custom check requested by the user
# For db_idx=49 (used in Fig 4), ω=1.0, new_ΔT=200.0.
# The prc2 value for db_idx=49 in Fig 3 was 0.41. A_ω_st_relation defaults to "NA" for prc2 if not given,
# which means it will use the default from the loaded parameter set for db_idx=49.
# Let's use the same PRC2 as Fig 3 for consistency if needed, otherwise let A_ω_st_relation use its default for db_idx 49
# Fig 3 used prc2 = 0.41 for db_idx = 49.
# Fig 4 uses db_idx = 49, but A_ω_st_relation is called without a prc2 override, so it uses the default for that db_idx from the database.
# For this custom check, we will also let it use the default prc2 for db_idx 49 by passing "NA" (which is default for A_ω_st_relation)
check_critical_amplitude_custom(db_idx=49, target_freq=1.0, new_ΔT=200.0, amplitude_range=0:5:300, prc2_val="NA")
check_critical_amplitude_custom(db_idx=49, target_freq=1.0, new_ΔT=100.0, amplitude_range=0:5:300, prc2_val="NA")

# Generate the new comparison figure for A-ω curves with different ΔT
generate_figure_compare_A_omega_for_different_deltaT()

# Print guidance message when the script is included or run without uncommenting a generation line
# Simplified condition to check if running as main script non-interactively
if !isinteractive() && abspath(PROGRAM_FILE) == @__FILE__
    println("\nScript finished setup. Check comments at the end for how to generate figures.")
    println("To generate figures:")
    println("  1. Uncomment a line like 'generate_figure3()' or 'run_all()' at the end of this script and run again.")
    println("  2. Or, run interactively: include(\"src/Notch_EMT_paper.jl\") then call functions like generate_figure3().")
end 
