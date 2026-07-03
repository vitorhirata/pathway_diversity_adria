include("src/common.jl")
using GeoMakie, GraphMakie, CairoMakie

# ── Constants ─────────────────────────────────────────────────────────────────

rcps = ["45", "70"]
dhw_scenarios = [5, 7, 11]
seed_years = 30
N_seed_total = 1e7
N_seed_weights = (
    N_seed_TA=0.15, N_seed_CA=0.5, N_seed_CNA=0.0, N_seed_SM=0.35, N_seed_LM=0.0
)

scenario_names = [
    :heat_stress, :geographic_spread, :connectivity,
    :functional_diversity, :balanced, :unguided
]

n_scenarios_per_dhw = 7  # 5 options + 1 unguided + 1 counterfactual
n_dhw = length(dhw_scenarios)
n_rcps = length(rcps)
n_options = length(scenario_names)  # 6 (excludes counterfactual)

# ── Domain setup (loaded once) ────────────────────────────────────────────────

dom = ADRIA.load_domain(pd_config["domain_path"], rcps[1];
    calib_params_fn=pd_config["coral_param_path"],
    timeframe=(2022, 2022 + seed_years + 2 + 10)
)
ms = ADRIA.model_spec(dom)

ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "FogCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "MCCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "Coral").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "GrowthAcceleration").fieldname)

ADRIA.fix_factor!(dom;
    wave_scenario=1,
    guided=1,
    min_iv_locations=200,
    plan_horizon=5.0,
    N_mc_settlers=0,
    fogging=0.0,
    SRM=0.0,
    seed_year_start=2,
    seed_years=seed_years,
    seed_deployment_freq=1,
    seed_strategy=1,
    seeding_devices_per_m2=5,
    a_adapt=5.0,
    a_adapt_ref=5,
    N_seed_TA=N_seed_total * N_seed_weights.N_seed_TA,
    N_seed_CA=N_seed_total * N_seed_weights.N_seed_CA,
    N_seed_CNA=N_seed_total * N_seed_weights.N_seed_CNA,
    N_seed_SM=N_seed_total * N_seed_weights.N_seed_SM,
    N_seed_LM=N_seed_total * N_seed_weights.N_seed_LM,
    depth_min=2.0,
    depth_offset=25.0
)

options = ADRIA.analysis.option_seed_preference(; include_weights=true)

# ── Build scenario table (direct column assignment, pawn pattern) ─────────────

n_scens_total = n_dhw * n_scenarios_per_dhw  # 3 × 7 = 21
scens = repeat(ADRIA.sample(dom, 2)[1:1, :], n_scens_total)

row = 1
for dhw in dhw_scenarios
    for option in eachrow(options)
        scens[row, :dhw_scenario] = dhw
        scens[row, :guided] = 1
        scens[row, :N_seed_TA] = N_seed_total * N_seed_weights.N_seed_TA
        scens[row, :N_seed_CA] = N_seed_total * N_seed_weights.N_seed_CA
        scens[row, :N_seed_CNA] = N_seed_total * N_seed_weights.N_seed_CNA
        scens[row, :N_seed_SM] = N_seed_total * N_seed_weights.N_seed_SM
        scens[row, :N_seed_LM] = N_seed_total * N_seed_weights.N_seed_LM
        scens[row, :seed_heat_stress] = option[2]
        scens[row, :seed_in_connectivity] = option[3]
        scens[row, :seed_out_connectivity] = option[4]
        scens[row, :seed_depth] = option[5]
        scens[row, :seed_coral_cover] = option[6]
        scens[row, :seed_cluster_diversity] = option[7]
        scens[row, :seed_geographic_separation] = option[8]
        scens[row, :seed_coral_diversity] = option[9]
        row += 1
    end

    # Unguided seeding
    scens[row, :dhw_scenario] = dhw
    scens[row, :guided] = 0
    scens[row, :N_seed_TA] = N_seed_total * N_seed_weights.N_seed_TA
    scens[row, :N_seed_CA] = N_seed_total * N_seed_weights.N_seed_CA
    scens[row, :N_seed_CNA] = N_seed_total * N_seed_weights.N_seed_CNA
    scens[row, :N_seed_SM] = N_seed_total * N_seed_weights.N_seed_SM
    scens[row, :N_seed_LM] = N_seed_total * N_seed_weights.N_seed_LM
    row += 1

    # Counterfactual (no seeding)
    scens[row, :dhw_scenario] = dhw
    scens[row, :guided] = 0
    scens[row, :N_seed_TA] = 0
    scens[row, :N_seed_CA] = 0
    scens[row, :N_seed_CNA] = 0
    scens[row, :N_seed_SM] = 0
    scens[row, :N_seed_LM] = 0
    row += 1
end

rs = ADRIA.run_scenarios(dom, scens, rcps)

# Load scenario
# path = "Outputs/"
# rs = ADRIA.load_results(path)

# ── Analysis parameters ───────────────────────────────────────────────────────

tail_fraction = 0.05
metric_labels = [
    "Years >20% coral cover\n(worst reefs)", "Years >20% coral cover\n(best reefs)",
    "Cumulative cover\n(worst reefs)", "Cumulative cover\n(best reefs)",
    "Cumulative evenness\n(worst reefs)", "Cumulative evenness\n(best reefs)"
]

"""
Select bottom and top `tail_fraction` of reefs by delta (opt − cf), then compute
(median(opt[group]) − median(cf[group])) / median(cf[group]) for each group.
Returns (bottom_ratio, top_ratio).
"""
function delta_tail_ratio(
    opt::AbstractVector, cf::AbstractVector; tail_fraction::Float64=0.05, norm::Float64 = 0.0
)
    delta = opt .- cf
    k = max(1, floor(Int, tail_fraction * length(delta)))
    order = sortperm(delta)
    bot = order[1:k]
    top = order[(end - k + 1):end]

    if iszero(norm)
        norm = mean(cf)
    end

    bot_ratio = mean(delta[bot]) / norm
    top_ratio = mean(delta[top]) / norm
    return bot_ratio, top_ratio
end

# ── Metrics over full timeframe (computed once for all scenarios) ──────────────

n_locs = size(rs.seed_log, :locations)
n_scens_rs = nrow(rs.inputs)  # n_scens_total × n_rcps after multi-RCP run

m_tac = Array(ADRIA.metrics.total_absolute_cover(rs)) .* 1e-6  # YAXArray (timesteps, locations, scenarios)
loc_hab_area_km2 = rs.loc_area .* rs.loc_max_coral_cover .* 1e-6
fd_data = ADRIA.metrics.coral_evenness(rs)               # YAXArray (timesteps, locations, scenarios)

n_yrs_above = ADRIA.ZeroDataCube((:locations, :scenarios), (n_locs, n_scens_rs); T=Int32)
for s in 1:n_scens_rs, l in 1:n_locs
    thr = 0.20 * loc_hab_area_km2[l]
    n_yrs_above.data[l, s] = Int32(count(m_tac[:, l, s] .>= thr))
end

cum_tac = dropdims(sum(m_tac; dims=1); dims=1)   # YAXArray (locations, scenarios)
cum_fd = dropdims(sum(fd_data; dims=:timesteps); dims=:timesteps)  # YAXArray (locations, scenarios)

# Reefs seeded in at least one intervention scenario (across all options, DHW, RCPs)
seed_start = Int(rs.inputs.seed_year_start[1])
n_seed_years = Int(rs.inputs.seed_years[1])
seed_ts = seed_start:(seed_start + n_seed_years - 1)
cf_idxs = n_scenarios_per_dhw:n_scenarios_per_dhw:n_scens_rs  # every 7th row is a counterfactual
intervention_idxs = setdiff(1:n_scens_rs, cf_idxs)
seed_per_reef_per_ts_scen = dropdims(
    sum(rs.seed_log[timesteps=seed_ts, scenarios=intervention_idxs]; dims=:coral_id);
    dims=:coral_id
)
never_seeded = vec(all(seed_per_reef_per_ts_scen.data .== 0; dims=(1, 3)))
active_mask = .!never_seeded

# ── CVaR aggregation ─────────────────────────────────────────────────────────

# all_results_aggreg_reef[rcp_i] is (n_options=6, n_dhw=3, n_metrics=6)
#   metrics 1-2: n_yrs_above lower/upper CVaR ratio vs counterfactual
#   metrics 3-4: cum_tac lower/upper CVaR ratio vs counterfactual
#   metrics 5-6: cum_fd lower/upper CVaR ratio vs counterfactual
# all_results_per_reef[rcp_i]    is (n_locs, n_options=6, n_dhw=3, n_metrics=3) — absolute deltas for boxplot
all_results_aggreg_reef = Vector{Array{Float64,3}}(undef, n_rcps)
all_results_per_reef = Vector{Array{Float64,4}}(undef, n_rcps)

for (rcp_i, rcp) in enumerate(rcps)
    rcp_idxs = findall(rs.inputs.RCP .== parse(Float64, rcp))
    results_aggreg_reef = Array{Float64}(undef, n_options, n_dhw, 6)
    results_per_reef = Array{Float64}(undef, n_locs, n_options, n_dhw, 3)

    for (d, _) in enumerate(dhw_scenarios)
        cf_s = rcp_idxs[(d - 1) * n_scenarios_per_dhw + n_scenarios_per_dhw]

        # Counterfactual absolute values over active reefs (same for all options in this block)
        cf_nyrs = Float64.(n_yrs_above.data[active_mask, cf_s])
        cf_tac  = cum_tac[active_mask, cf_s]
        cf_fd   = cum_fd.data[active_mask, cf_s]

        for (o, _) in enumerate(scenario_names)
            s = rcp_idxs[(d - 1) * n_scenarios_per_dhw + o]

            # Absolute values over active reefs
            opt_nyrs = Float64.(n_yrs_above.data[active_mask, s])
            opt_tac  = cum_tac[active_mask, s]
            opt_fd   = cum_fd.data[active_mask, s]

            # Absolute deltas over all reefs for boxplot
            results_per_reef[:, o, d, 1] = Float64.(n_yrs_above.data[:, s]) .- Float64.(n_yrs_above.data[:, cf_s])
            results_per_reef[:, o, d, 2] = cum_tac[:, s] .- cum_tac[:, cf_s]
            results_per_reef[:, o, d, 3] = cum_fd.data[:, s] .- cum_fd.data[:, cf_s]

            # Median ratio for bottom/top 10% reefs by delta
            results_aggreg_reef[o, d, 1], results_aggreg_reef[o, d, 2] =
                delta_tail_ratio(opt_nyrs, cf_nyrs; tail_fraction=tail_fraction, norm=float(seed_years))
            results_aggreg_reef[o, d, 3], results_aggreg_reef[o, d, 4] =
                delta_tail_ratio(opt_tac, cf_tac; tail_fraction=tail_fraction)
            results_aggreg_reef[o, d, 5], results_aggreg_reef[o, d, 6] =
                delta_tail_ratio(opt_fd, cf_fd; tail_fraction=tail_fraction)
        end
    end

    all_results_aggreg_reef[rcp_i] = results_aggreg_reef
    all_results_per_reef[rcp_i] = results_per_reef
end

# ── Visualisation ─────────────────────────────────────────────────────────────

option_colors = Makie.wong_colors()[1:n_options]
dhw_linestyles = [:solid, :dash, :dot]
dhw_markers = [:circle, :rect, :utriangle]
dhw_model_names = dom.dhw_scens.properties["model_names"][dhw_scenarios]

for (rcp_i, rcp) in enumerate(rcps)
    results = all_results_aggreg_reef[rcp_i]

    # ── Figure ───────────────────────────────────────────────────────────────

    fig = Figure(size=(1100, 480))
    ax = Axis(fig[1, 1];
        xticks=(1:6, metric_labels),
        xticklabelsize=11,
        ylabel="Relative performance against \ncounterfactual (no interv.)",
        title="Relative performance Δ/counterfactual - RCP $(rcp)"
    )

    hlines!(ax, [0]; color=:black, linewidth=2)
    vlines!(ax, [2.5, 4.5]; color=:gray70, linewidth=1, linestyle=:dash)

    for (o_i, _) in enumerate(scenario_names)
        color = option_colors[o_i]
        for (d_i, _) in enumerate(dhw_scenarios)
            y = results[o_i, d_i, :]
            lines!(ax, 1:6, y; color, linestyle=dhw_linestyles[d_i], linewidth=1.5)
            scatter!(ax, 1:6, y; color, marker=dhw_markers[d_i], markersize=10)
        end
    end

    Legend(fig[1, 2],
        [LineElement(; color=option_colors[i], linewidth=2) for i in 1:n_options],
        string.(scenario_names);
        title="Option", framevisible=false
    )

    Legend(fig[2, 2],
        [
            [LineElement(; linestyle=dhw_linestyles[i], color=:gray40, linewidth=2),
             MarkerElement(; marker=dhw_markers[i], color=:gray40)]
            for i in 1:n_dhw
        ],
        [dhw_model_names[i] for i in 1:n_dhw];
        title="Global Climate Model", framevisible=false
    )

    rowsize!(fig.layout, 2, Auto(0.4))

    save(
        joinpath(pd_config["plot_output_path"], "robustness_rcp$(rcp).png"),
        fig; px_per_unit=2
    )
    @info "Saved robustness_rcp$(rcp).png"
end

# ── Boxplot figures: per-reef distribution across DHW scenarios ───────────────

option_labels = string.(scenario_names)

for (rcp_i, rcp) in enumerate(rcps)
    per_reef = all_results_per_reef[rcp_i]  # (n_locs, n_options, n_dhw, 3)

    fig = Figure(size=(1200, 450))

    for (m_i, metric) in enumerate(metric_labels)
        ax = Axis(fig[1, m_i];
            xticks=(1:n_options, option_labels),
            xticklabelrotation=π / 4,
            ylabel=m_i == 1 ? "Relative difference vs counterfactual" : "",
            title=metric
        )
        hlines!(ax, [0]; color=:black, linewidth=2)

        for (o_i, _) in enumerate(scenario_names)
            # pool active reefs × all dhw_scenarios for this option
            y = vec(per_reef[active_mask, o_i, :, m_i])
            boxplot!(ax, fill(o_i, length(y)), y; color=option_colors[o_i])
        end
    end

    save(
        joinpath(pd_config["plot_output_path"], "robustness_boxplot_rcp$(rcp).png"),
        fig; px_per_unit=2
    )
    @info "Saved robustness_boxplot_rcp$(rcp).png"
end
