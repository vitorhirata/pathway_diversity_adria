#=
Convergence analysis. Try different sample sizes to analyse how much is necessary for pathway diversity results
to converge to its final value.
=#
include("src/common.jl")
using CairoMakie, Statistics, StatsBase

rs_path = ""
rs = ADRIA.load_results(rs_path)
scens = rs.inputs

# Select scenarios matching the first parameter combination in the result set
ref_N_seed_TA = scens.N_seed_TA[1]
ref_dhw_scenario = scens.dhw_scenario[1]
ref_min_iv_locs = scens.min_iv_locations[1]

condition = (
    scens.N_seed_TA .== ref_N_seed_TA .&&
    scens.dhw_scenario .== ref_dhw_scenario .&&
    scens.min_iv_locations .== ref_min_iv_locs
)
idx_scens_all = findall(condition)

# Reference pathway diversity using all matching scenarios (fraction = 1.0)
pd_ref = ADRIA.analysis.pathway_diversity(rs, idx_scens_all)

# Convergence loop: for each sample fraction, repeatedly sub-sample idx_scens_all
weights = 0.5:0.1:0.9
n_reps = 5

option_names = pd_ref.option_name
n_options = length(option_names)
n_rows = length(weights) * n_reps * n_options
results = DataFrame(;
    weight=zeros(Float64, n_rows),
    initial_option=fill("", n_rows),
    pathway_diversity=zeros(Float64, n_rows)
)

row_idx = 1
for (_, w) in enumerate(weights)
    @info "Running analysis for weight: $(w)"
    n_keep = ceil(Int, w * length(idx_scens_all))
    for _ in 1:n_reps
        idx_sample = sort(StatsBase.sample(idx_scens_all, n_keep; replace=false))
        removed_pathways = length(idx_scens_all) - n_keep
        pd_w = ADRIA.analysis.pathway_diversity(rs, idx_sample, removed_pathways)
        for pd_row in eachrow(pd_w)
            ref_val = pd_ref[pd_ref.option_name .== pd_row.option_name, :pathway_diversity][1]
            results[row_idx, :weight] = w
            results[row_idx, :initial_option] = String(pd_row.option_name)
            results[row_idx, :pathway_diversity] = pd_row.pathway_diversity - ref_val
            row_idx += 1
        end
    end
end

# Plot: one line per option, x = sample fraction, y = pd(w) - pd(1.0), band = min/max
fig = Figure(; size=(700, 450))
ax = Axis(fig[1, 1];
    xlabel="Sample fraction",
    ylabel="Pathway diversity − reference",
    xticks=collect(weights)
)

palette = Makie.current_default_theme().palette.color[]
for (i, opt) in enumerate(option_names)
    opt_rows = results[results.initial_option .== String(opt), :]
    xs = collect(weights)
    means = [median(opt_rows[opt_rows.weight .== w, :pathway_diversity]) for w in weights]
    lows = [
        quantile(opt_rows[opt_rows.weight .== w, :pathway_diversity], 0.1) for w in weights
    ]
    highs = [
        quantile(opt_rows[opt_rows.weight .== w, :pathway_diversity], 0.9) for w in weights
    ]
    col = palette[i]
    band!(ax, xs, lows, highs; color=(col, 0.2))
    lines!(ax, xs, means; color=col, label=string(opt))
end

hlines!(ax, [0.0]; color=:black, linestyle=:dash, linewidth=1)
axislegend(ax; position=:lb)
save(joinpath(pd_config["plot_output_path"], "pd_convergence.png"), fig; px_per_unit=2)
