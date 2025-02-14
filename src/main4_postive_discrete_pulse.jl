# This file is to explore the pulsatile effects on epigenetic switching.
## ====================== Loading packages and data library==========================
using Revise
using DifferentialEquations
using Plots;
gr(fontfamily="Helvetica");
# plotly()
using Catalyst
using Catalyst: parameters

using ProgressMeter
using DataFrames, CSV, Random, Distributions
using StatsPlots
using Latexify, Measures, FLoops, LaTeXStrings
includet("./Functions.jl")
db, p_names, initi_names, parameter_set, initial_condition, rand_idx_set, db_idx = loading_database()



## ============================= Import models =================================
signal_type = "bump"
model_bump = Import_model(; type=signal_type);
equations(model_bump);
signal_type = "pulsatile"
model_pulsatile = Import_model(; type=signal_type);
equations(model_pulsatile)
rn_latex, ode_latex = ode_latex_gen(model_pulsatile);

# latexify(model_pulsatile)
# Graph(model_pulsatile)
# ODE_equations = convert(ODESystem, model_pulsatile)
# @show latexify(ODE_equations
## ============================== single test run for the model 🍏 works =============================
# 🔴
db_idx = 592
freq = 0.0;
phase = 0;
amplitude = 220;
T_init = 1e-10; #! set this to be small
ΔT = 100;
tspan = (0.0, 350.0);
@show db[db_idx, p_names]
model = model_pulsatile

u0map, pmap, p, tspan, ts, cb, sol = single_solve(; db_idx=db_idx, freq=freq, phase=phase, amplitude=amplitude, T_init=T_init, ΔT=ΔT, tspan=tspan, phase_reset=true)
@show t_switching = find_ST(sol)

plt = plot(sol, idxs=[4, 5, 6, 9], lw=1.5, xlabel="Time", ylabel="Concentration", dpi=500)
osci_signal(t, A, w, ϕ) = A * (1 + sign(cos(w * t + ϕ)))
# osci_signal(t, A, w, ϕ) = A * (abs(cos(w * t + ϕ))) #🔴
tt = ts[1]:0.01:ts[2]
osci_signal.(tt, amplitude, freq, 0.0)
plot!(plt, tt, osci_signal.(tt, amplitude, freq, 0.0),
    label="Pulsatile Input", seriestype=:steppre, line=(:dot, 2), alpha=0.8,
    # ylims = [0, 700],
    fill=(0, 0.3, :darkgreen), color="black", dpi=300)





## ====== Changing the prc2 rate results in different swtiching times =========
prob_new = remake_prob(model_pulsatile, u0map, tspan, p; prc2=0.59)
@time sol = solve(prob_new, Rosenbrock23(), callback=cb, tstops=ts)
plt = plot(sol,
    idxs=[4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
    lw=1.5,
    xlabel="Time", ylabel="Concentration",
    title="PRC2 rate : 0.4",
    dpi=500)
anim_prc2 = anim_prc2_changing(0:0.1:0.8, tspan=[0, 350], u0=u0map)
## 


## ====== Single plot for 1 gene comparing pulsatile case vs bump case=======
T_init = 1e-10
# ======= pulsatile model
model = model_pulsatile
signal_type = "pulsatile"
plt_pulsatile = single_solve_plot(; model=model, db_idx=49, phase=0, freq=0.5, amplitude=65.0, T_init=T_init, ΔT=80, type=signal_type)
# ======== bump model
model = model_bump
signal_type = "bump"
plt_bump = single_solve_plot(; model=model, db_idx=49, phase=0, freq=0.5, amplitude=65.0, T_init=T_init, ΔT=80, type=signal_type)
plt_2model = plot(plt_pulsatile, plt_bump, layout=(2, 1))
##


# write a function to plot 1 gene comparing pulsatile case vs bump case
# # ! moved to single_gene_case.jl
# function single_solve_plot_pulsatile_bump(; db_idx=49, phase=0, freq=0.5, amplitude=65.0, T_init=0.01, ΔT=100)
#     plt_pulsatile = single_solve_plot(; model=model_pulsatile, db_idx=db_idx, phase=phase, freq=freq, amplitude=amplitude, T_init=T_init, ΔT=ΔT, type="pulsatile")
#     plt_bump = single_solve_plot(; model=model_bump, db_idx=db_idx, phase=phase, freq=freq, amplitude=amplitude, T_init=T_init, ΔT=ΔT, type="bump")
#     plot(plt_pulsatile, plt_bump, layout=(2, 1))
# end

# single_solve_plot_pulsatile_bump()


# #! Have to set T_init to 0.01 to avoid the discontinuity
# for ΔT ∈ 0.01:10:1000
#     plt = single_solve_plot_pulsatile_bump(ΔT=ΔT, T_init=0.01)
#     display(plt)
# end





## ====== single plot for 1 gene ID 592 ======= Dll4 vs Dll1 within the first the duartiona of pulse given T_init
T_init = 1e-10
# ======= sustained model
plt_sustained = single_solve_plot(; model=model_pulsatile, db_idx=592, phase=0, freq=0.0, amplitude=165.0, T_init=T_init, ΔT=50, type="sustained", phase_reset=true)
# ======= pulsatile model
plt_pulsatile = single_solve_plot(; model=model_pulsatile, db_idx=592, phase=0, freq=0.010, amplitude=165.0, T_init=T_init, ΔT=50, type="pulsatile", phase_reset=true)

plot(plt_sustained, plt_pulsatile, layout=(2, 1))
##


# * write function for ploting 1 gene ID default to 592 comparing sustained signal vs pulsatile signal.
function single_solve_plot_sustained_pulsatile(; db_idx=592, phase=0, freq=0.0, amplitude=165.0, T_init=0.01, ΔT=100)
    plt_sustained = single_solve_plot(; model=model_pulsatile, db_idx=db_idx, phase=0.0, freq=0.0, amplitude=amplitude, T_init=T_init, ΔT=ΔT, type="sustained", phase_reset=true)
    plt_pulsatile = single_solve_plot(; model=model_pulsatile, db_idx=db_idx, phase=phase, freq=freq, amplitude=amplitude, T_init=T_init, ΔT=ΔT, type="pulsatile", phase_reset=true)
    plot(plt_sustained, plt_pulsatile, layout=(2, 1))
end


# * showed some freq and amplitude allows pulsatile signal switch states after signal is off (between pulses)
for freq ∈ 0.01:0.05:1#, amplitude ∈ 165:10:300
    plt = single_solve_plot_sustained_pulsatile(ΔT=100, T_init=0.01, freq=freq, amplitude=300)#, amplitude = amplitude)
    display(plt)
    sleep(0.1)
end





## ======================== to find index with A = 64, Dll4 switch 
for id = 100:200
    find_id_Dll4_vs_Dll1(id, amplitude=300, prc2=0.64)
    sleep(0.1)
end
## ======================== to find index with A = 64, Dll4 switch





## * Exmaple 1 ======= Dll4 vs Dll1 for gene id:49 ====
    #! this gene displayed that the Dll1 signal can flip the histone state later then the first pulse comparing with Dll4 case.
single_gene_id = 49
id2_freq = 0.1535
amplitude1 = 62
amplitude2 = 62
T_init = 1e-10
ΔT = 100
prc2 = 0.64

plt_gene1_Dll4, plt_gene1_Dll1, plt_2genes_compare_id_49 =
    Two_Genes_TS_by_Prc2(;
        # model = model_bump,
        model=model_pulsatile,
        id1=single_gene_id, id2=single_gene_id,
        id2_freq=id2_freq, amplitude1=amplitude1,
        amplitude2=amplitude2, prc2=prc2,
        T_init=T_init, ΔT=ΔT, title_on=true, legend_title_on=false,
        vars_to_show=[5, 6, 9], #tspan_rt = 2, # 
        type="pulsatile",
        # type="bump",
        phase2_reset=true) #🔴 specify the type
plt_gene1_Dll4
plt_gene1_Dll1
plt_2genes_compare_id_49
# savefig(plt_2genes_compare_id_49,"./figures/APS_id_49_Dll4_Dll1_compare_bump.png")
##
# savepath = pathgen("pulsatile")
# savefig(plt_gene1_Dll4, savepath * "plt2_gene1_Dll4.png")
# savefig(plt_gene1_Dll1, savepath * "plt2_gene1_Dll1.png")




## * Example 2 ======= Dll4 vs Dll1 for gene id:592 ====
    #! this gene does not appear that activation time ST is different between Dll4 and Dll1. The activation is happened at the end of the first pulse
single_gene_id = 592
id2_freq = 0.09
amplitude1 = amplitude2 = 134
T_init = 1e-10
ΔT = 100
prc2 = 0.41

plt_gene2_Dll4, plt_gene2_Dll1, plt_2genes_compare_id_592 =
    Two_Genes_TS_by_Prc2(; model=model_pulsatile,
        id1=single_gene_id, id2=single_gene_id,
        id2_freq=id2_freq, amplitude1=amplitude1,
        amplitude2=amplitude2, prc2=prc2,
        T_init=T_init, ΔT=ΔT, title_on=true, legend_title_on=true,
        type="pulsatile",
        phase2_reset=true) #🔴 specify the type
plt_gene2_Dll4
plt_gene2_Dll1
plt_2genes_compare_id_592
##
# savepath = pathgen("bump")
# savefig(plt_gene2_Dll4, savepath * "plt_gene2_Dll4.png")
# savefig(plt_gene2_Dll1, savepath * "plt_gene2_Dll1.png")

## * Example 3 ======= Dll4 vs Dll1 for gene id:592 vs gene id:49 ========  
    #! this gene does not appear that activation time ST is different between Dll4 and Dll1. The activation is happened within the first pulse
gene_id_1 = 592
gene_id_2 = 49
id2_freq = 0.09
amplitude1 = amplitude2 = 104
T_init = 1e-10
ΔT = 100
prc2 = 0.11
plt_gene1_Dll4, plt_gene2_Dll1, plt_2genes_compare_id_592_49 =
    Two_Genes_TS_by_Prc2(; model=model_pulsatile,
        id1=gene_id_1, id2=gene_id_2,
        id2_freq=id2_freq, amplitude1=amplitude1,
        amplitude2=amplitude2, prc2=prc2,
        T_init=T_init, ΔT=ΔT, title_on=true, legend_title_on=true,
        type="pulsatile",
        phase2_reset=true)
plt_gene1_Dll4
plt_gene2_Dll1
plt_2genes_compare_id_592_49
## =============================================================================



## *======= Animation for gene ID:49 Dll4 vs Dll1, when varying PRC2 rates ===========  
gene_id = 49
id2_freq = 0.13
amplitude = 100
gr()
anim_comp = Two_genes_prc2_TS_animation(; prc2_range=0:0.02:1, model=model_pulsatile,
    id1=gene_id,
    id2=gene_id,
    id2_freq=id2_freq,
    amplitude1=amplitude,
    amplitude2=amplitude)
gif(anim_comp, fps=1)
##

## * gene ID: 49 as an example, visualization for the two genes TS by changing prc2 rate 

plt_Dll4, plt_Dll1, plt_2genes_compare = Two_Genes_TS_by_Prc2(; model=model_pulsatile,
    id1=49, id2=49, id2_freq=0.6,
    prc2=0.41,
    amplitude1=300, amplitude2=100,
    legend_title_on=false)
plt_2genes_compare
# plt_Dll4
# savefig(plt_Dll4, "./figures/Dll4_signal.png")
plt_Dll1
# savefig(plt_Dll1, "./figures/Dll1_signal_freq_1.4.png")

## =================================================================
#! --------------- case with multiple pulses -------------
T_init = 1e-10
# ======= pulsatile model
model = model_pulsatile
signal_type = "pulsatile"
plt_Dll1 = single_solve_plot(; model=model, db_idx=49, phase=0, freq=0.43, prc2=0.41, amplitude=50.0, T_init=T_init, ΔT=100, type=signal_type, T_final=1.5 * ΔT)
savefig(plt_pulsatile, "./figures/Dll1_multiple.png")

plt_Dll4 = single_solve_plot(; model=model, db_idx=49, phase=0, freq=0.0, prc2=0.41, amplitude=50.0, T_init=T_init, ΔT=100, type=signal_type, T_final=1.5 * ΔT)
savefig(plt_Dll4, "./figures/Dll4.png")
## =====================================================================





##

# savefig(plt_2genes_compare, "./figures/592_49_prc2=0.862.png")


## ------- loop over each parameter set to find large TS that is close to the ΔT[end]
# phase = 0
# freq = 0.0
# amplitude = 120
# T_init = 100
# ΔT = 100
# tspan = (0.0, 10^6)
# db_idx = 592
# for ΔT in 1:10^4:10^6
#     @show db_idx
#     @show ΔT
#     u0map, pmap, p, tspan, ts, cb, sol = single_solve(; db_idx = db_idx, freq = freq, phase = phase, amplitude = amplitude, T_init = T_init, ΔT = ΔT, tspan = tspan)
#     t_switching = switching_time(; sol = sol, pulse_period = T_init:0.1:T_init+ΔT, idx = [6, 9], return_plot = false)
#     @show t_switching
#     plt = plot(sol,
#         vars = [4, 5, 6, 9], # [MR,KDM5A,H4,H27,KDM6A]
#         lw = 2,
#         xlabel = "Time", ylabel = "Concentration",
#         foreground_color_legend = nothing,
#         title = "A : $amplitude, freq = $freq, switch time : $t_switching, ID : $db_idx",
#         dpi = 300)
#     display(plt)
# end

## -----🔶 animation of fixed amplitude and phase with increase frequency, and its switching time ----------------
# ------🔹 phase dependent case with phase ϕ = 0
# anim_freq_tswitch(; range=0:0.02:0.4, amplitude=400, db_idx=db_idx)




## 📗 ================  Calculate the A ω ϕ st relation for a single gene. specify the db_idx number for a gene.========================
db_idx = 49
# db_idx = 600 # 🍏 paper figure 6
df_save_path, figure_save_path = pathgen(db_idx=49, type="pulsatile")

# =============== Generate the the database for the A ω ϕ st relation =============================
amplitude_range = 0:1:300
freq_range = 0:0.02:2
df = A_ω_st_relation(; model=model_pulsatile,
    db_idx=db_idx,
    amplitude_range=amplitude_range,
    freq_range=freq_range,
    ΔT=100)
## ===========================================================================================


##======================  save df to df_save_path =================
filename = df_save_path * "freq :$freq_range" * "_|_" * "amplitude :$amplitude_range.csv"
# CSV.write(filename, df)
df = CSV.read(filename, DataFrame)
## ===========================================================================================

## ==== the relationship between driving frequency and switching time group by amplitude =========
plt_freq_vs_ST = df_freq_vs_ST_groupby_amp(df; amplitude_select=[], palette=cgrad([:goldenrod1, :dodgerblue4, :chartreuse3])) # show only 3 
#! =================== 6 representative amplitudes  ST - ω plots ------- Paper figure 6 (b) ==========================
pyplot()
plt_freq_vs_ST = df_freq_vs_ST_groupby_amp(df; amplitude_select=collect(50:50:300), figure_save_path=figure_save_path)

## ===========================================================================================


## ! ==== the relationship between driving amplitude and driving frequency A -ω , Paper figure 6 (a)========
# ---- first trim the database with each frequency and its associated minimum amplutide
df_min_amp = extract_min_amp(df)
# ---- plot the driving frequency vs driving minimum amplitude 🔴 Figure 6(a) -------
pyplot()
min_amp_vs_freq_plt = plot(df_min_amp.amp, df_min_amp.freq, 
    seriestype=:line,  
    linewidth=3,
    color=:red,
    legend=:topleft,
    label="Switching Boundary",
    xlabel=L"Driving Amplitude ($A$)", 
    ylabel=L"Driving Frequency ($\omega$)", 
    dpi=500)
savefig(min_amp_vs_freq_plt, figure_save_path * "freq_vs_amp.png")
## ===========================================================================================

## ======== Generating database for all gene id =====
function generate_database(; idx_range=nothing, amplitude_range=0:1:300, freq_range=0:0.02:2, prc2_range=0:0.05:1, ΔT=100)

    for db_idx ∈ idx_range
        df_save_path, _ = pathgen(db_idx=db_idx, type="pulsatile")
        println("Saving dataframe to $df_save_path")

        amplitude_range = amplitude_range
        freq_range = freq_range
        # df = A_ω_st_relation(; model=model_pulsatile,
        #     db_idx=db_idx,
        #     amplitude_range=amplitude_range,
        #     freq_range=freq_range,
        #     ΔT=ΔT)
        df = A_ω_st_relation_prc2_range(; model=model_pulsatile,
            db_idx=db_idx,
            amplitude_range=amplitude_range,
            freq_range=freq_range,
            prc2_range=prc2_range,
            ΔT=ΔT)
        savefile_name = df_save_path * "freq :$freq_range" * "_|_" * "amplitude :$amplitude_range.csv"
        println("Saving dataframe to $savefile_name")
        CSV.write(savefile_name, df)
    end
end




## ====================
#! working on this tonight
generate_database(idx_range=40:200,
    amplitude_range=50:50:300,
    freq_range=0:0.02:2,
    prc2_range=0:0.01:1,
    ΔT=100)
## ===========================================================================================







## =========== Generate the dataframe for A ω st prc2_range relation for a single gene. =========
db_idx = 592
amplitude_range = 0:100:300
freq_range = 0:0.1:2
df_592 = A_ω_st_relation_prc2_range(; model=model_pulsatile,
    db_idx=db_idx,
    amplitude_range=amplitude_range,
    freq_range=freq_range,
    prc2_range=0:0.1:1,
    ΔT=100)
## ===========================================================================================



##! ======================= Figure 7 A-w curve controled by prc2 rate ========================
db_idx = 592
amplitude_range = 0:1:500
freq_range = 0:0.05:1
df_592_3d = A_ω_st_relation_prc2_range(; model=model_pulsatile,
    db_idx=db_idx,
    amplitude_range=amplitude_range,
    freq_range=freq_range,
    prc2_range=0.1:0.1:1,
    ΔT=100)

# CSV.write("df_592_3d.csv", df_592_3d)
df_592_3d = CSV.read("df_592_3d.csv", DataFrame) # ! this is the dataframe for paper figure 7, can directly load it.
function plot_all_prc2(df::DataFrame, prc2_col::Symbol, x_col::Symbol, y_col::Symbol; kwargs...)
    # Use gr() backend instead of pyplot() since it's more stable
    pyplot()
    
    # group by prc2 column
    grouped_df = groupby(df, prc2_col)
    n_groups = length(grouped_df)
    
    # Create a color palette that transitions from cool to warm colors
    colors = cgrad(:viridis, n_groups, categorical=true)
    
    # create plot object
    plt = plot(xlabel="Driving Amplitude (A)",
        ylabel=L"Driving Frequency ($\omega$)", 
        dpi=500,
        legend=:topright,
        legendtitle="PRC2 Rate",
        legendfontsize=8,
        grid=false; kwargs...)
    
    # loop through each group and plot on same plot
    for (i, group) in enumerate(grouped_df)
        each_group = DataFrame(group)
        min_amp_df = extract_min_amp(each_group)
        
        # Plot only the line (no scatter points)
        plot!(plt, 
            min_amp_df[!, x_col], 
            min_amp_df[!, y_col], 
            linewidth=2,
            color=colors[i],
            label=string(round(each_group[1, prc2_col], digits=1)),
            legend_title="PRC2 Rate"
        )
    end
    
    return plt
end

A_w_prc2_curve = plot_all_prc2(df_592_3d, :prc2, :amp, :freq; legend=:outerright, foreground_color_legend=nothing, legendfontsize=7, legentitlefontsize=4, legendtitle="A-ω decision boundary\nregulated by PRC2 rate")
savefig(A_w_prc2_curve, "A_w_prc2_curve_id_592.png")


df_592_3d

# ===== added plotting function for 3d scatter plot for A-ω-stime groupby prc2 relation ======
# === New figure 🔴 ====
function plot_all_prc2_3d(df::DataFrame, prc2_col::Symbol, x_col::Symbol, y_col::Symbol, z_col::Symbol; kwargs...)
    # group by prc2 column
    grouped_df = groupby(df, prc2_col)
    # create plot object
    plt = plot(xlabel="Driving Amplitude (A)",
        ylabel="Driving Frequency (ω)", zlabel="Switching Time (ST)", dpi=500, titlefont=font(10, "Arial"))
    # loop through each group and plot on same plot with different colors
    for (i, group) in enumerate(grouped_df)
        # extract data for group
        each_group = DataFrame(group)
        # plot frequency vs amplitude vs switching time
    scatter!(plt, each_group[!, x_col], each_group[!, y_col], each_group[!, z_col], label=  string(each_group[1, prc2_col]); kwargs...)
    end
    return plt
end

plt_3d_prc2 = plot_all_prc2_3d(df_592_3d,
                 :prc2, :amp, :freq, :stime;
                m=1.5, guidefontsize=7, alpha = 0.6,
                legendtitle = " PRC2 rates",legend = :outerright)
savefig(plt_3d_prc2, "figures/df_592_3d_by_prc2.html")



## 🔴 ====== Generate the dataframe for the A ω st prc2 relation for a multi gene. ==========
df_49_592 = A_ω_st_prc2_multi_genes(gene_set_id=[49, 592],
    model=model_pulsatile,
    amplitude_range=50:50:300,
    freq_range=0:0.05:1,
    prc2_range=0:0.05:1,
    ΔT=100)
CSV.write("df_49_592.csv",df_49_592)
## ===========================================================================================




# show amplitude denpendent ST- ω relation for two genes.
#! this function not used in paper
function plot_ST_ω_by_amplitude_multi_genes(; data, layout=(3, 2), size=(1200, 800), legendfontsize=8, legendtitlefontsize=8, display_plt=false, prc2=0.3)
    fig_set = []
    for amplitude in unique(data.amp)
        fixed_A_2genes = filter(row -> row.amp == amplitude, data)
        fixed_A_prc2_2genes = filter(row -> row.prc2 == prc2, fixed_A_2genes)
        if !isempty(fixed_A_prc2_2genes)
            plt = @df fixed_A_prc2_2genes plot(
                :freq,
                :stime,
                group=:gene_id,
                palette=:tab10,
                m=(0.8, 1.5),
                # legend_title="Amplitude",
                legend_position=:outertopright,
                legendfontsize=legendfontsize,
                legendtitlefontsize=legendtitlefontsize,
                ylabel="Switching Time",
                xlabel="Switching Frequency",
                dpi=500,
                legend_title="Amplitude = $amplitude, \n PRC2  = $prc2"
            )
            push!(fig_set, plt)
        end
    end
    plt = plot(fig_set..., layout=layout, size=size)
    if display_plt
        display(plt)
    end

    return fig_set, plt
end

# plot_ST_ω_by_amplitude(data=df_49_592)
fig_set, plt = plot_ST_ω_by_amplitude_multi_genes(data=df_49_592)
plt







#! ================== The following figure is not used in paper yet, ST-ω curves with different A by different PRC2 values. ==========
function plot_freq_stime(df::DataFrame; selected_prc2::Vector{Float64}=Float64[], figsize=(800, 600))
    unique_prc2 = isempty(selected_prc2) ? unique(df.prc2) : selected_prc2
    unique_amp = unique(df.amp)
    unique_gene_id = unique(df.gene_id)
    n_prc2 = length(unique_prc2)
    layout = (ceil(Int, sqrt(n_prc2)), floor(Int, sqrt(n_prc2)))
    layout = layout[1] * layout[2] < n_prc2 ? (layout[1], layout[2] + 1) : layout
    plots = []
    gene_id_to_marker_shape = Dict(unique_gene_id[1] => :circle, unique_gene_id[2] => :cross)

    for prc2 in unique_prc2
        p = plot(xlabel="Frequency", ylabel="Switching Time", title="PRC2 = $prc2", legend=:topleft)

        for gene_id in unique_gene_id
            marker_shape = gene_id_to_marker_shape[gene_id]

            for amp in unique_amp
                sub_df = df[(df.prc2.==prc2).&(df.gene_id.==gene_id).&(df.amp.==amp), :]
                plot!(p, sub_df.freq, sub_df.stime, label="A = $amp", markershape=marker_shape, seriestype=:line, legend=:outertopright)
            end
        end

        push!(plots, p)
    end

    return plot(plots..., layout=layout, size=figsize)
end

plot_freq_stime(df_49_592, selected_prc2=[0.2, 0.4, 0.6, 0.8], figsize=(1200, 800))





# ******************⬇ Used section for multi genes ⬇ ********************
## ==== #! ============== Figure 9 , 6 plots for amplitude-dependent ST_ω relation for 2 genes ==============
df_49_592# ! new dataframe
fixed_amp = 100
gene2_Dll4_df = filter(row -> row.freq == 0.1 && row.amp == fixed_amp && row.gene_id == 592, df_49_592)
gene1_Dll1_df = filter(row -> row.freq == 0.1 && row.amp == fixed_amp && row.gene_id == 49, df_49_592)

plt1 = @df gene1_Dll1_df plot(:prc2, :stime,
    # palette=:RdYlBu_6,
    palette=:tab10,
    # m=(0.8, 2),
    legend_title="Genes",
    legend_position=:topright,
    linewidth=2,
    xlabel="PRC2 rate",
    ylabel="Switching Time (ST)",
    dpi=500,
    foreground_color_legend=nothing,
    title="Amplitude = $fixed_amp",
    label="Gene 1 (Pusatile frequency = 0.9)",
)
plt2 = @df gene2_Dll4_df plot!(plt1, :prc2, :stime,
    palette=:RdYlBu_6,
    # m=(0.8, 2),
    legend_title="Genes",
    legend_position=:topright,
    linewidth=2,
    xlabel="PRC2 rate",
    ylabel="Switching Time (ST)",
    dpi=500,
    foreground_color_legend=nothing,
    title="Amplitude = $fixed_amp",
    label="Gene 2 (Sustained)"
)
# ylims!(plt2, (100, 150))
# savefig(plt2, "plt2_ST_ω_2gene_prc2_increase_fixed_amp_400_freq_0.9.png")
##




## ------- plot ω vs. ST for various prc2 with amplitude  = 300
df_49_592 #! new df
df_49 = filter(row -> row.gene_id == 49, df_49_592)
df_49_fixed_amp = filter(row -> row.amp == 300, df_49)

plt_id49_amp_200_prc2_set = @df df_49_fixed_amp plot(
    :freq,
    :stime,
    group=:prc2,
    palette=palette(:inferno,22),  # Changed to a more diverse palette to highlight middle range values more distinctly
    m=(1, 2.5),
    legend_position=:outertopright,
    legendfontsize=5,
    legendtitlefontsize=5,
    ylabel="Switching Time",
    xlabel="Driving Frequency",
    dpi=500,
    legend_title="Amplitude = 200, \n PRC2 rate"
)
savefig(plt_id49_amp_200_prc2_set, "plt_id49_amp_200_prc2_set.png")
## 



## ------- plot ω vs. ST for various prc2 with amplitude  = 100
df_49 = filter(row -> row.gene_id == 49, df_49_592)
df_49_fixed_amp = filter(row -> row.amp == 300, df_49)
plt_id49_amp_100_prc2_set = @df df_49_fixed_amp plot(
    :freq,
    :stime,
    group=:prc2,
    palette=:RdBu_9,
    m=(1, 2.5),
    lw=2,
    # legend_title="Amplitude",
    legend_position=:outertopright,
    legendfontsize=5,
    legendtitlefontsize=5,
    ylabel="Switching Time",
    xlabel="Driving Frequency",
    dpi=500,
    legend_title="Amplitude = 100, \n PRC2 rate"
)
savefig(plt_id49_amp_100_prc2_set, "plt_id49_amp_100_prc2_set.png")
## 







## =================================================================
#! loading two dataframes
# Example usage
path = joinpath(dirname(@__DIR__), "Data/regular")
prefix = "db_idx:"
idx_range = 40:60  # Replace with the desired range of indices
csv_filename = "freq :0.0:0.02:2.0_|_amplitude :50:50:300.csv"
df_set, combined_df = load_csv_files(path, prefix, idx_range, csv_filename)
combined_df
## =================================================================



## ================================================= Select two genes for comparision ===========================================
df_2genes = filter_by_gene_id(combined_df, rand(40:60, 2))
df_2genes = filter_by_gene_id(combined_df, [6, 8])# !working example 
df_2genes = filter_by_gene_id(combined_df, [43, 52])# !working example 
#! working cases: [ 43, 52]
# [ 48, 55] : A = 50 cross, A = 300, no cross
# =================================================================================

for Dll1_freq in 0.1:0.1:0.9
    plt_prc2_ST = plot_prc2_vs_ST_amp_dependent(df_2genes=df_2genes,
        fixed_amp_set=[50, 100, 150, 200, 250, 300],
        freq_selection=[Dll1_freq, 0],
        size=(800, 800),
        fontsize=font(8),
        xlims=(0, 1))
    display(plt_prc2_ST)
    sleep(0.5)
end


## ============================display and saving figure 9 ======================================
#! choose two amplidues and plot for Figure 9
fixed_amp_set = df_2genes.amp |> unique
plot_Dll1vs_Dll4_prc2_ST_by_amplitude(; df_2genes=df_2genes, fixed_amp_set=fixed_amp_set, freq_selection=[0.5, 0], layout=(3, 2), size=(1000, 800), fontsize=font(8), save=true, filename="Dll1vs_Dll4_prc2_ST_by_amplitude_6amplitudes")
## ==================================================================

# ******************⬆ Used section for multi genes ⬆ ********************

##  plot the dynamics plot of gene id 43 with amplitude 300 and frequency 0.5
T_init = 1e-10
# ======= pulsatile model
model = model_pulsatile
signal_type = "pulsatile"

# !before the critial prc2 rate which is  around 0.35
plt_Dll1_gene1 = single_solve_plot(; model=model, db_idx=49, phase=0, freq=0.5, prc2=0.2, amplitude=300.0, T_init=T_init, ΔT=20, type=signal_type, T_final=1.5 * ΔT)

plt_Dll1_gene2 = single_solve_plot(; model=model, db_idx=52, phase=0, freq=0., prc2=0.2, amplitude=300.0, T_init=T_init, ΔT=20, type=signal_type, T_final=1.5 * ΔT)

# !after the critial prc2 rate which is  around 0.35
plt_Dll1_gene1_2 = single_solve_plot(; model=model, db_idx=49, phase=0, freq=0.5, prc2=0.5, amplitude=300.0, T_init=T_init, ΔT=20, type=signal_type, T_final=1.5 * ΔT)

plt_Dll1_gene2_2 = single_solve_plot(; model=model, db_idx=52, phase=0, freq=0., prc2=0.5, amplitude=300.0, T_init=T_init, ΔT=20, type=signal_type, T_final=1.5 * ΔT)


## =====
## 🔴 ====== Generate the dataframe for the A ω st prc2 relation for a multi gene. ==========
df_43_52 = A_ω_st_prc2_multi_genes(gene_set_id=[43, 52],
    model=model_pulsatile,
    amplitude_range=50:50:300,
    freq_range=0:0.1:1,
    prc2_range=0:0.05:1,
    ΔT=100)
# CSV.write("df_43_52.csv", df_43_52)
# !== load the above generated dataframe df_43_52 ========
df_43_52 = CSV.read("df_43_52.csv", DataFrame)


# --- filter gene 43, with amplitude 300, and frequency 0.5 and prc2 0.2
gene1_ST_b = filter(row -> row.gene_id == 43 && row.amp == 300 && row.freq == 0.5 && row.prc2 == 0.2, df_43_52)

# ----filter gene 52 with amplitude 300, and frequency 0 and prc2 0.2
gene2_ST_b = filter(row -> row.gene_id == 52 && row.amp == 300 && row.freq == 0.0 && row.prc2 == 0.2, df_43_52)

# ----filter gene 43 with amplitude 300, and frequency 0.5 and prc2 0.4
gene1_ST_a = filter(row -> row.gene_id == 43 && row.amp == 300 && row.freq == 0.5 && row.prc2 == 0.4, df_43_52)

# ----filter gene 52 with amplitude 300, and frequency 0 and prc2 0.4
gene2_ST_a = filter(row -> row.gene_id == 52 && row.amp == 300 && row.freq == 0.0 && row.prc2 == 0.4, df_43_52)


## ----------------
gene_id_1 = 52
gene_id_2 = 43
id2_freq = 0.5
amplitude1 = amplitude2 = 150
T_init = 1e-10
ΔT = 30
prc2 = 0.2
plt_gene1_Dll4_before, plt_gene2_Dll1_before, plt_2genes_compare_id_52_43_before =
    Two_Genes_TS_by_Prc2(; model=model_pulsatile,
        id1=gene_id_1, id2=gene_id_2,
        id2_freq=id2_freq, amplitude1=amplitude1,
        amplitude2=amplitude2, prc2=prc2,
        T_init=T_init, ΔT=ΔT, title_on=true, legend_title_on=true,
        type="pulsatile",
        phase2_reset=true)
plt_2genes_compare_id_52_43_before
savefig(plt_2genes_compare_id_52_43_before, "plt_2genes_compare_id_52_43_before.png")
## ----
prc2 = 0.4
plt_gene1_Dll4_after, plt_gene2_Dll1_after, plt_2genes_compare_id_52_43_after =
    Two_Genes_TS_by_Prc2(; model=model_pulsatile,
        id1=gene_id_1, id2=gene_id_2,
        id2_freq=id2_freq, amplitude1=amplitude1,
        amplitude2=amplitude2, prc2=prc2,
        T_init=T_init, ΔT=ΔT, title_on=true, legend_title_on=true,
        type="pulsatile",
        phase2_reset=true)
plt_2genes_compare_id_52_43_after
savefig(plt_2genes_compare_id_52_43_after, "plt_2genes_compare_id_52_43_after.png")



## ---- twicking gene1 with freq = 0.5----- 
prc2 = 0.2
gene1_Dll4, gene1_Dll1, plt_test =
    Two_Genes_TS_by_Prc2(; model=model_pulsatile,
        id1=gene_id_1, id2=gene_id_1,
        id2_freq=0.5, amplitude1=amplitude1,
        amplitude2=amplitude2, prc2=prc2,
        T_init=T_init, ΔT=9, title_on=true, legend_title_on=true,
        type="pulsatile",
        phase2_reset=true)
plt_test
plot(gene1_Dll1, legend_title="Gene 1")
savefig(plt_test, "plt_2genes_compare_id_52_43_1stgene_freq_0.5.png")

plt_gene1_Dll4_before, plt_gene2_Dll1_before, plt_2genes_compare_id_52_43_before =
    Two_Genes_TS_by_Prc2(; model=model_pulsatile,
        id1=gene_id_1, id2=gene_id_2,
        id2_freq=id2_freq, amplitude1=amplitude1,
        amplitude2=amplitude2, prc2=prc2,
        T_init=T_init, ΔT=9, title_on=true, legend_title_on=true,
        type="pulsatile",
        phase2_reset=true)
plt_2genes_compare_id_52_43_before
twogenes_pulsatile = plot(plot(gene1_Dll1, legend_title="Gene 1"), plt_gene2_Dll1_before, layout=(2, 1))
savefig(twogenes_pulsatile, "plt_2genes_compare_id_52_43_1stgene_freq_0.5.png")

## --- twicking gene2 with freq = 0.1 -------
prc2 = 0.4
_,_, plt_test =
    Two_Genes_TS_by_Prc2(; model=model_pulsatile,
        id1=gene_id_1, id2=gene_id_2,
        id2_freq=0.1, amplitude1=amplitude1,
        amplitude2=amplitude2, prc2=prc2,
        T_init=T_init, ΔT=9, title_on=true, legend_title_on=true,
        type="pulsatile",
        phase2_reset=true)
plt_test
savefig(plt_test, "plt_2genes_compare_id_52_43_2ndgene_freq_0.1.png")
## ---

