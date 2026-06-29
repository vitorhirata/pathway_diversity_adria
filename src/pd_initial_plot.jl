#=
Script showing results to support pathway diversity simulation.
Runs a simplified ADRIA simulation where all parameters are fixed except the seeding option,
over a 10-year horizon with a single decision step where the option can change. Plots the
switching probability matrix, all six of its subcomponents and the probability distribution
for that decision step.
=#

include("src/common.jl")
using CairoMakie, GeoMakie, GraphMakie

# ----------------------------------------------------------
# Constants
RCP = "45"        # RCP 26, 45, 70
dhw_scenario = 11
N_seed = 1e7        # total seeds; split across taxa via N_seed_weights below
min_iv_locations = 200
seed_year_start = 2
seed_years = 10
pd_frequency = 5
N_seed_weights = (
    N_seed_TA=0.15, N_seed_CA=0.5, N_seed_CNA=0.0, N_seed_SM=0.35, N_seed_LM=0.0
)

# ----------------------------------------------------------
# Build & run the simplified simulation
dom = ADRIA.load_domain(
    pd_config["domain_path"], RCP;
    calib_params_fn=pd_config["coral_param_path"],
    # timeframe: seed_years + 2 (start seeding). Final timestep = 2022 + seed_years + 2
    timeframe=(2022, 2022 + seed_years + 2)
)
ms = ADRIA.model_spec(dom)

ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "FogCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "MCCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "Coral").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "GrowthAcceleration").fieldname)

ADRIA.fix_factor!(dom;
    #Seeding params
    seed_year_start=seed_year_start, # Start as soon as possible
    seed_years=seed_years,           # Pathway diversity analysis time
    seed_deployment_freq=1,          # Seed every year
    seeding_devices_per_m2=5,
    seed_strategy=1,                 # Periodic deployment
    a_adapt=5.0,
    a_adapt_ref=5,
    # Interventions params
    plan_horizon=5.0,
    guided=1,                        # CoCoSo
    #Decision params
    depth_min=2.0,
    depth_offset=25.0,
    #Other interventions (no fogging, shading or moving corals)
    fogging=0.0,
    SRM=0.0,
    N_mc_settlers=0,
    #Seeding amounts
    N_seed_TA=N_seed * N_seed_weights.N_seed_TA,
    N_seed_CA=N_seed * N_seed_weights.N_seed_CA,
    N_seed_CNA=N_seed * N_seed_weights.N_seed_CNA,
    N_seed_SM=N_seed * N_seed_weights.N_seed_SM,
    N_seed_LM=N_seed * N_seed_weights.N_seed_LM,
    #Environmental params
    dhw_scenario=dhw_scenario,
    min_iv_locations=min_iv_locations,
    wave_scenario=1
)

# Generate scenarios
scens = ADRIA.sample_options(dom, pd_frequency)

# Run scenarios
rs = ADRIA.run_scenarios(dom, scens, RCP)

# ----------------------------------------------------------
# Load scenarios
# path = ""
# rs = ADRIA.load_results(path)

# ----------------------------------------------------------
# Options definition plot
figures_path = joinpath(pd_config["plot_output_path"], "pd_figures")
mkpath(figures_path)
options_iw = ADRIA.analysis.option_seed_preference(; include_weights=true)
fig = ADRIA.viz.mcda_options(options_iw)
save(joinpath(pd_config["plot_output_path"], "pd_figures", "option_weight_matrix.png"), fig)

# ----------------------------------------------------------
# Single-decision-step matrices
options = ADRIA.analysis.option_seed_preference()  # 5 options, defines row/col order
option_names = options.option_name
n_opt = length(option_names)
opt_idx = Dict(o => i for (i, o) in enumerate(option_names))

idx_scens = collect(1:nrow(scens))
max_time = size(rs.seed_log, :timesteps)
decoded_ts = Dict(
    i => ADRIA.analysis.decode_option_ts(
        rs.inputs.option_ts[i], seed_year_start, seed_years, pd_frequency, max_time
    )
    for i in idx_scens
)
tstep = seed_year_start + pd_frequency  # single decision step
mcda_method = ADRIA.mcda_methods()[Int64(rs.inputs.guided[1])]

# Outcome-perf caches over window [tstep, t_end]
rel_cover = Array(rs.outcomes[:relative_cover][scenarios=idx_scens])
fd_data = Array(rs.outcomes[:coral_evenness][scenarios=idx_scens])
k_area = ADRIA.loc_k_area(rs)
scen_to_idx = Dict(s => i for (i, s) in enumerate(idx_scens))
t_end = min(tstep + pd_frequency - 1, max_time)
# Per-option performance: one representative scenario per option (holding that option at
# tstep), with per-location cumulative metrics (metric × location YAXArray).
option_perf = ADRIA.analysis._compute_option_perf(
    rel_cover, fd_data, k_area, idx_scens, scen_to_idx, decoded_ts, tstep, t_end
)

dms = ADRIA.readcubedata(
    rs.decision_matrix_log[timesteps=ADRIA.At([tstep]), scenarios=idx_scens]
)
ports = ADRIA.analysis._ports()
# Valid locations are those represented in the logged decision matrix
valid_locs = collect(1:size(rs.loc_data, 1))

prob_matrix = zeros(n_opt, n_opt)
cum_rel_tac_diff = zeros(n_opt, n_opt)
cum_fd_diff = zeros(n_opt, n_opt)
distance_port = zeros(n_opt, n_opt)
dispersion = zeros(n_opt, n_opt)
similarity = zeros(n_opt, n_opt)

for i in idx_scens
    src = decoded_ts[i][tstep - 1]
    dst = decoded_ts[i][tstep]
    r = opt_idx[src]
    c = opt_idx[dst]
    dm = dms[timesteps=ADRIA.At(tstep), scenarios=ADRIA.At(i)]

    # combined (normalized) switching probability to dst -- matches the model
    prob_matrix[r, c] = ADRIA.analysis.switching_probability(
        src, dm, rs.loc_data, mcda_method, min_iv_locations, dst, option_perf
    )

    # selected locations for src/dst from THIS scenario's decision matrix
    src_locs = ADRIA.select_locations(
        options[options.option_name .== src, :preference][1], dm, mcda_method,
        valid_locs, min_iv_locations
    )
    dst_locs = ADRIA.select_locations(
        options[options.option_name .== dst, :preference][1], dm, mcda_method,
        valid_locs, min_iv_locations
    )
    src_df = ADRIA.analysis._filter_locations(rs.loc_data, src_locs)
    dst_df = ADRIA.analysis._filter_locations(rs.loc_data, dst_locs)
    common = intersect(src_df.UNIQUE_ID, dst_df.UNIQUE_ID)
    u_src = src_df[src_df.UNIQUE_ID .∉ [common], :]
    u_dst = dst_df[dst_df.UNIQUE_ID .∉ [common], :]

    # 3 spatial subcomponents
    distance_port[r, c] = ADRIA.analysis.distance_port_score(u_dst, u_src, ports)
    dispersion[r, c] = ADRIA.analysis.dispersion_score(dst_df, src_df)
    similarity[r, c] = ADRIA.analysis.option_similarity(u_dst, u_src)

    # 2 outcome subcomponents (mirror switching_probability internals: per-location
    # paired difference between candidate and past option representatives)
    past_performance = option_perf[src]
    candidate_performance = option_perf[dst]

    cum_rel_tac_diff[r, c] = ADRIA.analysis.two_sided_cvar(
        candidate_performance[metric=ADRIA.At(:rel_tac)] ./ past_performance[metric=ADRIA.At(:rel_tac)] .- 1;
        σ=ADRIA.analysis._σ_rel_tac
    )
    cum_fd_diff[r, c] = ADRIA.analysis.two_sided_cvar(
        candidate_performance[metric=ADRIA.At(:fd)] ./ past_performance[metric=ADRIA.At(:fd)] .- 1;
        σ=ADRIA.analysis._σ_fd
    )
end

# ----------------------------------------------------------
# Switching probability matrix plot
option_labels = string.(option_names)
fig = Figure()
ax = Axis(
    fig[1, 1];
    xlabel="To",
    ylabel="From",
    xticks=(1:n_opt, option_labels),
    yticks=(1:n_opt, option_labels),
    xticklabelrotation=2.0 / π,
    yreversed=true
)
limits = (10.0, 30.0)
heatmap!(ax, transpose(100 * prob_matrix); colorrange=limits)
Colorbar(fig[1, 2]; limits=limits)
text!(ax,
    string.(round.(vec(100 * prob_matrix); digits=1));
    position=[Point2f(x, y) for x in 1:n_opt for y in 1:n_opt],
    align=(:center, :center),
    color=:white,
    fontsize=14
)
save(joinpath(pd_config["plot_output_path"], "pd_figures", "probability_matrix.png"), fig)

# ----------------------------------------------------------
# Subcomponents grid: switching probability + 6 subcomponents
function _heatmap_panel!(ax, mat, colorrange=(0.0, 1.0), digits=2)
    heatmap!(ax, transpose(mat); colorrange=colorrange)
    return text!(ax,
        string.(round.(vec(mat); digits=digits));
        position=[Point2f(x, y) for x in 1:n_opt for y in 1:n_opt],
        align=(:center, :center),
        color=:white,
        fontsize=12
    )
end

function _panel_axis(fig, pos, title; ylabel="")
    return Axis(
        fig[pos...];
        title=title,
        xlabel="To",
        ylabel=ylabel,
        xticks=(1:n_opt, option_labels),
        yticks=(1:n_opt, option_labels),
        xticklabelrotation=2.0 / π,
        yreversed=true
    )
end

fig = Figure(; size=(3 * 520, 2 * 450))
# Row 1
ax1 = _panel_axis(fig, (1, 1), "Switching probability (%)"; ylabel="From")
heatmap!(ax1, transpose(prob_matrix); colorrange=(0.0, 1.0))
text!(ax1,
    string.(round.(vec(100 * prob_matrix); digits=1));
    position=[Point2f(x, y) for x in 1:n_opt for y in 1:n_opt],
    align=(:center, :center),
    color=:white,
    fontsize=12
)
_heatmap_panel!(_panel_axis(fig, (1, 2), "Cumulative relative TAC diff"), cum_rel_tac_diff)
_heatmap_panel!(_panel_axis(fig, (1, 3), "Cumulative functional diversity diff"), cum_fd_diff)
# Row 2
_heatmap_panel!(_panel_axis(fig, (2, 1), "Distance to port"; ylabel="From"), distance_port)
_heatmap_panel!(_panel_axis(fig, (2, 2), "Dispersion"), dispersion)
_heatmap_panel!(_panel_axis(fig, (2, 3), "Option similarity"), similarity)
Colorbar(fig[1:2, 4]; limits=(0.0, 1.0))
save(joinpath(pd_config["plot_output_path"], "pd_figures", "option_scores_matrix.png"), fig)

# ----------------------------------------------------------
# 1-D plot: mean over rows for each destination option (col)
fig = Figure()
ax = Axis(
    fig[1, 1];
    xlabel="Option",
    ylabel="Mean score",
    xticks=(1:n_opt, option_labels),
    xticklabelrotation=π / 4,
    yticks=(0.0:0.2:1.0),
    limits=(nothing, (0.0, 1.0))
)
scatterlines!(ax, 1:n_opt, vec(mean(cum_rel_tac_diff; dims=1)); label="Cum. rel. TAC diff")
scatterlines!(ax, 1:n_opt, vec(mean(cum_fd_diff; dims=1)); label="Cum. FD diff")
scatterlines!(ax, 1:n_opt, vec(mean(distance_port; dims=1)); label="Distance to port")
scatterlines!(ax, 1:n_opt, vec(mean(dispersion; dims=1)); label="Dispersion")
scatterlines!(ax, 1:n_opt, vec(mean(similarity; dims=1)); label="Option similarity")
scatterlines!(
    ax, 1:n_opt, vec(mean(prob_matrix; dims=1)); color=:black, label="Switching prob"
)
fig[1, 2] = Legend(fig, ax; framevisible=false)
save(joinpath(pd_config["plot_output_path"], "pd_figures", "option_scores.png"), fig)

# ----------------------------------------------------------
# Probability distribution histogram
n_bins = 50
all_probs = ADRIA.analysis.pathway_diversity(rs, idx_scens; scenario_probabilities=true).probability
# Shared bin edges so every panel uses the same x-scale and aligned bins
edges = range(minimum(minimum.(all_probs)), maximum(maximum.(all_probs)); length=n_bins + 1)

fig = Figure(; size=(1600, 400))
axes = Axis[]
for (idx, option) in enumerate(option_names)
    ax = Axis(fig[1, idx]; xlabel=string(option))
    hist!(ax, all_probs[idx]; bins=edges)
    push!(axes, ax)
end
linkaxes!(axes...)  # same x and y scale across all panels
save(
    joinpath(pd_config["plot_output_path"], "pd_figures", "probabilities_histogram.png"),
    fig
)
