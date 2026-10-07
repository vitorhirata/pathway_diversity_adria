using CSV, DataFrames, CairoMakie

#=
#### Setup for collecting data ####
1. Run pathway_diversity_analysis.jl code defining params for the correct set of parameters to be plotted
2. Load a result_set: rs = ADRIA.load_results(path)
3. Remove Threads.@threads and add a log on ADRIA.jl pathway_diversity.jl. function delta_tail_ratio before return
    _raw_file = σ == _σ_rel_tac ? "./Outputs/raw_rel_tac.csv" : "./Outputs/raw_fd.csv"
    open(_raw_file, "a") do io
        println(io, "$(cvar_lower),$(cvar_upper),$(raw)")
    end
=#

param = params[1]
rcp="45"
condition =
    rs.inputs.N_seed_CA .== param[1] * N_seed_weights.N_seed_CA .&&
    rs.inputs.N_seed_TA .== param[1] * N_seed_weights.N_seed_TA .&&
    rs.inputs.min_iv_locations .== param[2] .&&
    rs.inputs.dhw_scenario .== param[3] .&&
    rs.inputs.RCP .== parse(Float64, rcp)
idx_scens = findall(condition)
scenario_result = ADRIA.analysis.pathway_diversity(rs, idx_scens)

# Rename files to add number of seeds and run again for another set of parameters

param = params[4]
rcp="45"
condition =
    rs.inputs.N_seed_CA .== param[1] * N_seed_weights.N_seed_CA .&&
    rs.inputs.N_seed_TA .== param[1] * N_seed_weights.N_seed_TA .&&
    rs.inputs.min_iv_locations .== param[2] .&&
    rs.inputs.dhw_scenario .== param[3] .&&
    rs.inputs.RCP .== parse(Float64, rcp)
idx_scens = findall(condition)
scenario_result = ADRIA.analysis.pathway_diversity(rs, idx_scens)

# Rename files again, load data and plot results

raw_rel_tac_low_seed = CSV.read("./Outputs/raw_rel_tac_seed_1e6.csv", DataFrame; header=false, types=Float64)
raw_rel_tac_high_seed = CSV.read("./Outputs/raw_rel_tac_seed_1e8.csv", DataFrame; header=false, types=Float64)

raw_fd_low_seed = CSV.read("./Outputs/raw_fd_seed_1e6.csv", DataFrame; header=false, types=Float64)
raw_fd_high_seed = CSV.read("./Outputs/raw_fd_seed_1e8.csv", DataFrame; header=false, types=Float64)

raw_rel_tac = (raw_rel_tac_low_seed, raw_rel_tac_high_seed)
raw_fd = (raw_fd_low_seed, raw_fd_high_seed)

titles = ("N_corals=1M, N_reefs=100", "N_corals=100M, N_reefs=500")

raw = -0.002:0.0001:0.002
scales = [0.0005, 0.001, 0.002]
f(x, scale) = 0.5 * (1 + tanh(x / scale))

fig = Figure(; size=(1200, 1200))

# Add a bold panel label (e.g. "(a)") to the top-left of an axis
panel_label!(ax, label) = text!(
    ax, 0.02, 0.98; text=label, font=:bold, align=(:left, :top), space=:relative
)

# Row 1: cumulative cover histograms
for (i, data) in enumerate(raw_rel_tac)
    ax = Axis(fig[1, i]; title=titles[i], xlabel="cum. cover best + worst", ylabel=(i == 1 ? "count" : ""))
    hist!(ax, data.Column3; bins=40, color=(:steelblue, 0.7))
    panel_label!(ax, "($(('a':'e')[i]))")
end

# Row 2: cumulative evenness histograms
for (i, data) in enumerate(raw_fd)
    ax = Axis(fig[2, i]; title=titles[i], xlabel="cum. evenness best + worst", ylabel=(i == 1 ? "count" : ""))
    hist!(ax, data.Column3; bins=40, color=(:darkorange, 0.7))
    panel_label!(ax, "($(('c':'e')[i]))")
end

# Row 3: benefit function, spanning both columns
ax = Axis(fig[3, 1:2]; xlabel="performance metric best + worst", ylabel="f(x)")
for scale in scales
    lines!(ax, raw, [f(x, scale) for x in raw]; label="scale = $scale")
end
axislegend(ax; position=:rb)
panel_label!(ax, "(e)")

save("./Outputs/benefit_distribution.png", fig)
