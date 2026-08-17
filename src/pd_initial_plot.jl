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
dhw_scenario = 7
N_seed = 1e7        # total seeds; split across taxa via N_seed_weights below
min_iv_locations = 300
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
fix_common_parameters!(dom)

ADRIA.fix_factor!(dom;
    #Seeding params
    seed_years=seed_years,           # Pathway diversity analysis time
    seeding_devices_per_m2=5,
    a_adapt=5.0,
    #Seeding amounts
    N_seed_TA=N_seed * N_seed_weights.N_seed_TA,
    N_seed_CA=N_seed * N_seed_weights.N_seed_CA,
    N_seed_CNA=N_seed * N_seed_weights.N_seed_CNA,
    N_seed_SM=N_seed * N_seed_weights.N_seed_SM,
    N_seed_LM=N_seed * N_seed_weights.N_seed_LM,
    #Environmental params
    dhw_scenario=dhw_scenario,
    min_iv_locations=min_iv_locations
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

# Outcome-perf caches over window [tstep, t_end].
# Use per-location total absolute cover (km²) rather than relative cover for the outcome metric.
abs_cover_data = Array(ADRIA.metrics.total_absolute_cover(rs)[scenarios=idx_scens]) .* 1e-6
fd_data = Array(rs.outcomes[:coral_evenness][scenarios=idx_scens])
scen_to_idx = Dict(s => i for (i, s) in enumerate(idx_scens))
t_end = min(tstep + pd_frequency - 1, max_time)
decision_times = collect(seed_year_start:pd_frequency:(seed_year_start + seed_years - 1))
perf_options = vcat(option_names, :nothing)
# Per-prefix performance: keyed by (past_prefix..., candidate_option) tuples. At the first
# switching decision (tstep), keys are 2-tuples (src_option, candidate_option).
perf_by_prefix = ADRIA.analysis._prefix_option_perf(
    abs_cover_data, fd_data, idx_scens, scen_to_idx, decoded_ts, decision_times, tstep, t_end
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
similarity = zeros(n_opt, n_opt)

for i in idx_scens
    src = decoded_ts[i][tstep - 1]
    dst = decoded_ts[i][tstep]
    dst == :nothing && continue  # skip do-nothing scenarios (not plotted)
    r = opt_idx[src]
    c = opt_idx[dst]
    dm = dms[timesteps=ADRIA.At(tstep), scenarios=ADRIA.At(i)]

    # Build per-src option_perf dict (mirrors pathway_diversity internals): keyed by
    # candidate option symbol, performance conditioned on the same past prefix (src).
    option_perf::Dict{Symbol, YAXArray} = Dict(
        o => perf_by_prefix[(src, o)]
        for o in perf_options
        if haskey(perf_by_prefix, (src, o))
    )

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

    # 2 spatial subcomponents
    distance_port[r, c] = ADRIA.analysis.distance_port_score(u_dst, u_src, ports)
    similarity[r, c] = ADRIA.analysis.option_similarity(u_dst, u_src)

    # 2 outcome subcomponents (mirror switching_probability internals): keeping the same
    # option scores 0.7 (inertia bias); switching compares candidate against :nothing counterfactual.
    if src == dst
        cum_rel_tac_diff[r, c] = 0.7
        cum_fd_diff[r, c] = 0.7
    elseif haskey(option_perf, :nothing) && haskey(option_perf, dst)
        counter = option_perf[:nothing]
        cand_perf = option_perf[dst]
        cum_rel_tac_diff[r, c] = ADRIA.analysis.delta_tail_ratio(
            cand_perf[metric=ADRIA.At(:rel_tac)] ./ counter[metric=ADRIA.At(:rel_tac)] .- 1;
            σ=ADRIA.analysis._σ_rel_tac
        )
        cum_fd_diff[r, c] = ADRIA.analysis.delta_tail_ratio(
            cand_perf[metric=ADRIA.At(:fd)] ./ counter[metric=ADRIA.At(:fd)] .- 1;
            σ=ADRIA.analysis._σ_fd
        )
    end
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

function _panel_axis(fig, pos, title; ylabel="", xlabel="", kwargs...)
    return Axis(
        fig[pos...];
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        xticks=(1:n_opt, option_labels),
        yticks=(1:n_opt, option_labels),
        xticklabelrotation=2.0 / π,
        yreversed=true,
        titlesize=22,
        xlabelsize=22,
        ylabelsize=22,
        xticklabelsize=18,
        yticklabelsize=18,
        kwargs...
    )
end

fig = Figure()
panel_px = 350   # fixed square size shared by every heatmap panel
cbar_width = 25

# Left: switching probability heatmap, same size as the others and vertically centred
ax1 = _panel_axis(
    fig, (1:2, 2), "Switching probability (%)";
    ylabel="From", width=panel_px, height=panel_px, valign=:center
)
heatmap!(ax1, 100 .* transpose(prob_matrix); colorrange=(10.0, 32.0))
text!(ax1,
    string.(round.(vec(100 * prob_matrix); digits=1));
    position=[Point2f(x, y) for x in 1:n_opt for y in 1:n_opt],
    align=(:center, :center),
    color=:white,
    fontsize=12
)
# Left colorbar for the switching probability, matched to the panel height
Colorbar(
    fig[1:2, 1]; limits=(10.0, 32.0), flipaxis=false, width=cbar_width, height=panel_px,
    valign=:center, label="Switching probability", labelsize=22, ticklabelsize=18
)

# Right: 2x2 grid of subcomponents, aligned with each other and matched to the panel size
_heatmap_panel!(
    _panel_axis(
        fig, (1, 3), "Cumulative absolute cover";
        ylabel="From", width=panel_px, height=panel_px, xticklabelsvisible=false
    ),
    cum_rel_tac_diff
)
_heatmap_panel!(
    _panel_axis(
        fig, (1, 4), "Cumulative functional diversity";
        width=panel_px, height=panel_px, xticklabelsvisible=false, yticklabelsvisible=false
    ),
    cum_fd_diff
)
_heatmap_panel!(
    _panel_axis(fig, (2, 3), "Distance to port"; ylabel="From", xlabel="To", width=panel_px, height=panel_px, ),
    distance_port
)
_heatmap_panel!(
    _panel_axis(
        fig, (2, 4), "Option similarity"; xlabel="To", width=panel_px, height=panel_px, yticklabelsvisible=false
    ),
    similarity
)
# Right colorbar, matched in size to the left one
Colorbar(
    fig[1:2, 5]; limits=(0.0, 1.0), width=cbar_width, height=panel_px, valign=:center,
    ticklabelsize=18, label="Scores used in probability", labelsize=22
)
resize_to_layout!(fig)
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
all_probs = filter(!iszero, all_probs)
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
