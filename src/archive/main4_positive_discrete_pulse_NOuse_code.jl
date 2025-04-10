## now test for how prc2 changes will affect 2 genes ST_ω relation in a amplitude denpendent manner
# ---- generate dataframe for 2 genes (49, 592) as prc2 rate increases ------
#! =================== ⬇ this function is abandon, no ϕ now. ⬇=================
function Gen_df_stack_prc2_increase(; model=model, amplitude_range=0:50:300, freq_range=0:0.01:2, gene_set_id=[49, 592], phase_sample_size=10, prc2_range=0:0.1:1, mute_parameter_disp=true)
    df_stack_prc2_set = []
    @showprogress for prc2 in prc2_range
        df_stack_prc2_i = A_ω_ϕ_st_relation_multi_gene(amplitude_range=amplitude_range, freq_range=freq_range, gene_set_id=gene_set_id, phase_sample_size=phase_sample_size, prc2=prc2, mute_parameter_disp=mute_parameter_disp)
        if isempty(df_stack_prc2_i)
            println("No data for prc2 = ", prc2)
            continue
        else
            df_stack_prc2_i[!, :prc2] .= prc2
            @show df_stack_prc2_i
        end
        push!(df_stack_prc2_set, df_stack_prc2_i)
    end
    return vcat(df_stack_prc2_set...)
end



# ======= generate data for two genes with various prc2=======
df_stack_prc2_set = Gen_df_stack_prc2_increase(; model=model_pulsatile, amplitude_range=400, freq_range=0:0.1:2, gene_set_id=[49, 592], phase_sample_size=15, prc2_range=0.3:0.001:0.4, mute_parameter_disp=true)
CSV.write("df_stack_prc2_set_amplitude_range=400, freq_range=0:0.1:2, gene_set_id=[49, 592], phase_sample_size=15, prc2_range=0.3:0.001:0.4.csv", df_stack_prc2_set)
df_stack_prc2_set = CSV.File("old/data/df_stack_prc2_set_amplitude_range=300, freq_range=0:0.1:2, gene_set_id=[49, 592], phase_sample_size=15, prc2_range=0.3:0.001:0.4.csv") |> DataFrame
show(df_stack_prc2_set, allrows=true)
unique(df_stack_prc2_set.gene_id)
df_stack_prc2_set_new_gene_name = df_stack_prc2_set
# filter!(row -> row.amp == 300, df_stack_prc2_set_new_gene_name)

# generate a more detailed data for finner prc2 values.
# df_stack_prc2_finner_set = Gen_df_stack_prc2_increase(; amplitude_range=50:50:300, freq_range=0:0.05:2, gene_set_id=[49, 592], phase_sample_size=15, prc2_range=0.3:0.1:0.6, mute_parameter_disp=true)
# df_stack_prc2_finner_set
# CSV.write("df_stack_prc2_finner_set_amplitude_range=50:50:300, freq_range=0:0.05:2, gene_set_id=[49, 592], phase_sample_size=15, prc2_range=0.3:0.01:0.6.csv", df_stack_prc2_finner_set)
#! ================ ⬆ this function is abandon, no ϕ now.⬆ =================



##! ===================== ⬇ Abandoned Code ⬇ =====================
# =============== save and laod thie dataframe to a .csv file
using CSV
# CSV.write("df_stack_prc2_set_amplitude_range=0:50:400, freq_range=0:0.01:2, gene_set_id=[49, 592], phase_sample_size=10, prc2_range=0:0.1:1.csv", df_stack_prc2_set)
df_stack_prc2_set = CSV.File("../old/data/df_stack_prc2_set_amplitude_range=0:50:400, freq_range=0:0.01:2, gene_set_id=[49, 592], phase_sample_size=10, prc2_range=0:0.1:1.csv") |> DataFrame

# ============ the loaded dataframe has gene_id with 49 and 592, I just want to map this column to "gene_1" and "gene_2"
df_stack_prc2_set_new_gene_name = df_stack_prc2_set
idx_gne_dict = Dict(49 => "Gene 1", 592 => "Gene 2")
df_stack_prc2_set_new_gene_name = transform(df_stack_prc2_set_new_gene_name, :gene_id => ByRow(x -> idx_gne_dict[x]) => :gene_id)
unique(df_stack_prc2_set_new_gene_name.gene_id)
gdf_stack_prc2 = groupby(df_stack_prc2_set_new_gene_name, :prc2)

plot_ST_ω_by_amplitude(; data=gdf_stack_prc2[3], layout=(3, 3), size=(2200, 1600), display_plt=true)




# ----compare ST ω relation only at fixed amplitude
# check prc2 values 
function select_df_prc2_fixed_amp(; data=gdf_stack_prc2, prc2_key=3, amp_select=300)
    df_prc2 = data[keys(data)[prc2_key]]
    df_prc2_fixed_amp = filter(row -> row.amp == amp_select, df_prc2)
    println("returning prc2 = ", keys(data)[prc2_key], " and amp = ", amp_select)
    return df_prc2_fixed_amp
end





keys(gdf_stack_prc2)[9]
fixed_amp = 300
df_stack_prc2_fixed_amp = select_df_prc2_fixed_amp(prc2_key=3, amp_select=fixed_amp)
fig_set, plt = plot_ST_ω_by_amplitude(data=df_stack_prc2_fixed_amp, layout=(1, 1), size=(900, 600))
plt


## ====== fixed amplitude  change prc2 for two genes =================
plt_set = []
fixed_amp = 300
for prc2_idx in 1:length(keys(gdf_stack_prc2))
    @show prc2_idx
    df_stack_prc2_fixed_amp = select_df_prc2_fixed_amp(prc2_key=prc2_idx, amp_select=fixed_amp)
    fig_set, plt = plot_ST_ω_by_amplitude(data=df_stack_prc2_fixed_amp,
        layout=(1, 1), size=(900, 600),
        legendfontsize=5, legendtitlefontsize=5,
        display_plt=false, prc2=keys(gdf_stack_prc2)[prc2_idx][1])
    # title!(plt, "PRC2 = $(keys(gdf_stack_prc2)[prc2_idx][1])", titlefont=font(8))
    xlims!(plt, (0, 2))
    ylims!(plt, (100, 150))
    push!(plt_set, plt)
end

plt_ST_ω_2gene_prc2_increase_fixed_amp = plot(plt_set[1:4]..., layout=(2, 2), size=(1600, 1000), legendfontsize=12, legendtitlefontsize=12, m=(3, 0.6), legend=:topright)
# savefig(plt_ST_ω_2gene_prc2_increase_fixed_amp, "plt_ST_ω_2gene_prc2_increase_fixed_amp.png")




# ========== make an animation for prc2 0.3 ~ 0.4 at amp = 300 for two genes (49,592) ==================
plt_set = []
fixed_amp = 300
anim_prc2 = @animate for prc2_idx in 1:length(keys(gdf_stack_prc2))
    @show prc2_idx
    df_stack_prc2_fixed_amp = select_df_prc2_fixed_amp(prc2_key=prc2_idx, amp_select=fixed_amp)
    fig_set, plt = plot_ST_ω_by_amplitude(data=df_stack_prc2_fixed_amp,
        layout=(1, 1), size=(900, 600),
        legendfontsize=5, legendtitlefontsize=5,
        display_plt=false, prc2=keys(gdf_stack_prc2)[prc2_idx][1])
    # title!(plt, "PRC2 = $(keys(gdf_stack_prc2)[prc2_idx][1])", titlefont=font(8))
    xlims!(plt, (0, 2))
    ylims!(plt, (100, 200))
    push!(plt_set, plt)
end

gif(anim_prc2, "prc2_0.3~0.4_amp_200_id(49,592).gif", fps=15)


fixed_amp = 200
gene2_sustained_df = filter(row -> row.freq == 0 && row.gene_id == "Gene 2", df_stack_prc2_set_new_gene_name)
gene1_pusatile_df = filter(row -> row.freq == 0.9 && row.gene_id == "Gene 1", df_stack_prc2_set_new_gene_name)

##! ============================⬆ Abandoned Code⬆ =====================


## ========= To generate a plot of ω vs. ST for a range of PRC2 rate and a range of amplitude ==========
# ------- generate dataset
# df_stack_prc2_set_gene_49 = Gen_df_stack_prc2_increase(; model=model_pulsatile, amplitude_range=100:200:300, freq_range=0:0.1:2, gene_set_id=[49], phase_sample_size=15, prc2_range=0.1:0.05:0.6, mute_parameter_disp=true)
df_stack_prc2_set_gene_49 = Gen_df_stack_prc2_increase(; model=model_pulsatile, amplitude_range=100:200:300, freq_range=0:0.1:2, gene_set_id=[49], phase_sample_size=5, prc2_range=0.41, mute_parameter_disp=true)
# CSV.write("df_stack_prc2_set_gene_49.csv", df_stack_prc2_set_gene_49)



# ## * ===## =============== A-w-ϕ curve if ϕ is not controllable================
# switch_amplitude = []
# switch_frequency = []
# T_init = 1e-10
# tspan = (0.0, 600.0)
# ΔT = 100.0
# @showprogress for freq in 0.0:0.1:1#exp10.(-4:0.05:1)
#     for amplitude = 0:10:2000#15:0.1:20#14:0.01:18#10:0.05:30
#         switch_set = []
#         for phase = 0:2π
#             # frequency modulation
#             # p1 = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, freq, 0.0]
#             # p = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, freq, 0.0, phase]
#             # u0 = [6.0, 0.0, 6.0, 40.0, 500.0, 0.6876, 0.678, 500.0, 50.6344, 1.0, 2.0]
#             # u0 = [6.59, 0.0, 6.59, 43.41, 61.47, 0.02, 0.54, 48.17, 49.43, 1.09, 1.09]
#             # @show parameter_set[db_idx, :]
#             # @show initial_condition[db_idx, :]

#             p = vcat([collect(parameter_set[db_idx, :]), freq, 0.0, phase]...)
#             pmap = parameters(model) .=> p
#             u0 = collect(initial_condition[db_idx, :])
#             u0map = species(model) .=> u0
#             ts, cb = make_cb([T_init, T_init + ΔT], 13, amplitude)
#             prob1 = ODEProblem(model, u0map, tspan, pmap)
#             sol1 = solve(prob1, Rosenbrock23(), callback=cb, tstops=ts)

#             # plt1 = plot(sol1, vars = [5, 6, 9, 10], lw = 1.5, title = "Amp: $amplitude, Frequency: $freq")
#             # display(plt1)

#             check = check_switching(sol1, ts, tspan)
#             append!(switch_set, check)
#             # if check == -1 # if there exist a phase
#             #     append!(switch_amplitude, amplitude)
#             #     append!(switch_frequency, freq)
#             #     # append!(switch_phase, phase)
#             # end
#             if check_switching(sol1, ts, tspan) == -1
#                 break
#             end
#         end
#         # @show switch_set
#         if -1 in switch_set
#             append!(switch_amplitude, amplitude)
#             append!(switch_frequency, freq)
#         end
#         if -1 in switch_set
#             break
#         end
#     end
# end