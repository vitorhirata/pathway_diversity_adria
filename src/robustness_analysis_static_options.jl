include("src/common.jl")
using GeoMakie, GraphMakie, CairoMakie

# ── Constants ─────────────────────────────────────────────────────────────────

rcps = ["45", "70"]
dhw_scenarios = [2, 7, 10]
seed_years = 20
N_seed_totals = [1e6, 1e8]
min_iv_locations_list = [100, 200, 500]
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
    # Baseline only; every scenario row overrides min_iv_locations explicitly below.
    min_iv_locations=min_iv_locations_list[1],
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
    # Baseline only; every scenario row overrides N_seed_* explicitly below.
    N_seed_TA=N_seed_totals[1] * N_seed_weights.N_seed_TA,
    N_seed_CA=N_seed_totals[1] * N_seed_weights.N_seed_CA,
    N_seed_CNA=N_seed_totals[1] * N_seed_weights.N_seed_CNA,
    N_seed_SM=N_seed_totals[1] * N_seed_weights.N_seed_SM,
    N_seed_LM=N_seed_totals[1] * N_seed_weights.N_seed_LM,
    depth_min=2.0,
    depth_offset=25.0
)

options = ADRIA.analysis.option_seed_preference(; include_weights=true)

# ── Build scenario table (direct column assignment, pawn pattern) ─────────────

n_nseed = length(N_seed_totals)
n_minloc = length(min_iv_locations_list)
# Row layout (outer → inner): N_seed_total → min_iv_locations → dhw →
# [5 options, unguided, counterfactual]. Every n_scenarios_per_dhw-th row is a counterfactual.
n_scens_total = n_nseed * n_minloc * n_dhw * n_scenarios_per_dhw  # 2 × 3 × 3 × 7 = 126
scens = repeat(ADRIA.sample(dom, 2)[1:1, :], n_scens_total)
# Identifier carried through to rs.inputs: -1 = counterfactual, 0 = unguided, 1:5 = guided
# option in the same order as `options` / scenario_names[1:5].
scens[!, :option] = zeros(Int, n_scens_total)

row = 1
for N_seed_total in N_seed_totals, min_iv in min_iv_locations_list, dhw in dhw_scenarios
    for (opt_i, option) in enumerate(eachrow(options))
        scens[row, :dhw_scenario] = dhw
        scens[row, :guided] = 1
        scens[row, :option] = opt_i
        scens[row, :min_iv_locations] = min_iv
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
    scens[row, :option] = 0
    scens[row, :min_iv_locations] = min_iv
    scens[row, :N_seed_TA] = N_seed_total * N_seed_weights.N_seed_TA
    scens[row, :N_seed_CA] = N_seed_total * N_seed_weights.N_seed_CA
    scens[row, :N_seed_CNA] = N_seed_total * N_seed_weights.N_seed_CNA
    scens[row, :N_seed_SM] = N_seed_total * N_seed_weights.N_seed_SM
    scens[row, :N_seed_LM] = N_seed_total * N_seed_weights.N_seed_LM
    row += 1

    # Counterfactual (no seeding) — identical across N_seed/min_iv, duplicated to keep 7-row blocks
    scens[row, :dhw_scenario] = dhw
    scens[row, :guided] = 0
    scens[row, :option] = -1
    scens[row, :min_iv_locations] = min_iv
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

tail_fraction = 0.03

"""
Select bottom and top `tail_fraction` of reefs by delta (opt − cf), then normalize each tail
by its own counterfactual mean. Returns `(bottom_ratio, top_ratio, raw)` where
`raw = top_ratio − |bottom_ratio|` (the net variation, mirroring
`two_sided_cvar`'s `cvar_upper − |cvar_lower|`).

Note: this is a *ratio-of-means* (`mean(delta[tail]) / mean(cf[tail])`), whereas
`two_sided_cvar` is a *mean-of-ratios* (`mean(opt_l/cf_l − 1)`). They agree when cf is roughly
flat across the tail and diverge where a tail reef has near-zero cf. `norm` is a zero-guard
fallback used when a tail's counterfactual mean is ~0 (relevant for `n_yrs_above`, where many
counterfactual reefs never exceed the 20% threshold).
"""
function delta_tail_ratio(
    opt::AbstractVector, cf::AbstractVector; tail_fraction::Float64=0.03, norm::Float64 = 0.0
)
    delta = opt .- cf
    k = max(1, floor(Int, tail_fraction * length(delta)))
    order = sortperm(delta)
    bot = order[1:k]
    top = order[(end - k + 1):end]

    bot_ratio = mean(delta[bot]) / (iszero(norm) ? mean(cf) : norm)
    top_ratio = mean(delta[top]) / (iszero(norm) ? mean(cf) : norm)

    raw = top_ratio - abs(bot_ratio)
    return bot_ratio, top_ratio, raw
end

# ── Metrics: base quantities computed once (horizon-independent) ──────────────

n_locs = size(rs.seed_log, :locations)
n_scens_rs = nrow(rs.inputs)  # n_scens_total × n_rcps after multi-RCP run

# `readcubedata` materialises the lazy YAXArrays into memory while keeping them YAXArrays
# (dims preserved), so the named slicing below stays cheap. (`read` would strip dims to an Array.)
m_tac = ADRIA.readcubedata(ADRIA.metrics.total_absolute_cover(rs)) .* 1e-6  # YAXArray (timesteps, locations, scenarios), km²
fd_arr = ADRIA.readcubedata(ADRIA.metrics.coral_evenness(rs))               # YAXArray (timesteps, locations, scenarios)
loc_hab_area_km2 = rs.loc_area .* rs.loc_max_coral_cover .* 1e-6

# Seed window start (measurement horizons are anchored here).
seed_start = Int(rs.inputs.seed_year_start[1])

# ── Visualisation styling (shared across horizons) ────────────────────────────

option_colors = Makie.wong_colors()[1:n_options]
dhw_linestyles = [:solid, :dash, :dot]
dhw_markers = [:circle, :rect, :utriangle]
dhw_model_names = dom.dhw_scens.properties["model_names"][dhw_scenarios]

# Aggregate figure: 3 metrics × {worst, best, net} = 9 columns
aggreg_metric_labels = [
    "Years >20%\n(worst)", "Years >20%\n(best)", "Years >20%\n(net)",
    "Cum. cover\n(worst)", "Cum. cover\n(best)", "Cum. cover\n(net)",
    "Cum. evenness\n(worst)", "Cum. evenness\n(best)", "Cum. evenness\n(net)"
]
boxplot_metric_labels = ["Years >20% coral cover", "Cumulative cover", "Cumulative evenness"]
option_labels = string.(scenario_names)

# ── CVaR aggregation + plotting (per N_seed, RCP, horizon) ────────────────────

# Flat metric coordinates for the returned YAXArrays.
aggreg_metric_syms = [
    :nyrs_worst, :nyrs_best, :nyrs_net,
    :tac_worst, :tac_best, :tac_net,
    :fd_worst, :fd_best, :fd_net
]
base_metric_syms = [:nyrs, :tac, :fd]

"""
    cvar_aggregation(rs, m_tac, fd_arr, loc_hab_area_km2,
                     N_seed_total, min_iv, rcp, horizon; tail_fraction=0.03)

Filter the static-option scenarios for a given `N_seed_total`, `min_iv` (min_iv_locations),
`rcp` and measurement `horizon` (years from seed-start), window the metrics, and
CVaR-aggregate each option vs its counterfactual. Returns `(aggreg, per_reef)`:
- `aggreg`   : YAXArray `(option=6, dhw_scenario=3, metric=9)` — (worst, best, net) tail ratios
  vs counterfactual (over all reefs) for n_yrs_above / cum cover / cum evenness.
- `per_reef` : YAXArray `(option=6, dhw_scenario=3, metric=3, location=n_locs)` — absolute
  Δ vs counterfactual per location, for the 3 base metrics.
"""
function cvar_aggregation(
    rs, m_tac, fd_arr, loc_hab_area_km2,
    N_seed_total, min_iv, rcp, horizon; tail_fraction=0.03
)
    @assert seed_start + horizon - 1 <= size(m_tac, :timesteps) "Run too short for horizon $(horizon)."
    window = seed_start:(seed_start + horizon - 1)

    rcp_val = parse(Float64, rcp)
    # Common filter for this parameter set; the counterfactual carries no N_seed so it is
    # matched by `param_mask` alone, while seeded rows additionally match the seed budget.
    param_mask = (rs.inputs.RCP .== rcp_val) .& (rs.inputs.min_iv_locations .== min_iv)
    seeded_mask = param_mask .& (rs.inputs.N_seed_CA .== N_seed_total * N_seed_weights.N_seed_CA)

    # Windowed per-location metrics for a single scenario column (named-dim slicing)
    win_nyrs(s) = Float64[
        count(m_tac[timesteps=window, locations=l, scenarios=s].data .>= 0.20 * loc_hab_area_km2[l])
        for l in 1:n_locs
    ]
    win_tac(s) = vec(sum(m_tac[timesteps=window, scenarios=s].data; dims=1))
    win_fd(s)  = vec(sum(fd_arr[timesteps=window, scenarios=s].data; dims=1))

    aggreg = ADRIA.DataCube(
        zeros(n_options, n_dhw, 9);
        option=scenario_names, dhw_scenario=dhw_scenarios, metric=aggreg_metric_syms
    )
    per_reef = ADRIA.DataCube(
        zeros(n_options, n_dhw, 3, n_locs);
        option=scenario_names, dhw_scenario=dhw_scenarios,
        metric=base_metric_syms, location=1:n_locs
    )

    for (d, dhw) in enumerate(dhw_scenarios)
        dhw_mask = rs.inputs.dhw_scenario .== dhw
        # Counterfactual (option = -1) carries no seed budget, so match on param_mask only.
        cf_s = findfirst(param_mask .& dhw_mask .& (rs.inputs.option .== -1))

        cf_nyrs = win_nyrs(cf_s)
        cf_tac  = win_tac(cf_s)
        cf_fd   = win_fd(cf_s)

        for (o, _) in enumerate(scenario_names)
            # option id: 1:5 for guided options, 0 for the unguided (6th) column.
            opt_id = o < n_options ? o : 0
            s = findfirst(seeded_mask .& dhw_mask .& (rs.inputs.option .== opt_id))
            opt_nyrs = win_nyrs(s)
            opt_tac  = win_tac(s)
            opt_fd   = win_fd(s)

            # Per-reef absolute deltas over all locations (boxplot)
            per_reef[option=o, dhw_scenario=d, metric=ADRIA.At(:nyrs)] .= opt_nyrs .- cf_nyrs
            per_reef[option=o, dhw_scenario=d, metric=ADRIA.At(:tac)]  .= opt_tac  .- cf_tac
            per_reef[option=o, dhw_scenario=d, metric=ADRIA.At(:fd)]   .= opt_fd   .- cf_fd

            # (worst, best, net) tail ratios vs counterfactual over all reefs. n_yrs uses
            # horizon length as the zero-guard fallback (many cf reefs never exceed threshold).
            aggreg[option=o, dhw_scenario=d, metric=ADRIA.At([:nyrs_worst, :nyrs_best, :nyrs_net])] .= delta_tail_ratio(
                opt_nyrs, cf_nyrs; tail_fraction=tail_fraction, norm=float(horizon)
            )
            aggreg[option=o, dhw_scenario=d, metric=ADRIA.At([:tac_worst, :tac_best, :tac_net])] .= delta_tail_ratio(
                opt_tac, cf_tac; tail_fraction=tail_fraction
            )
            aggreg[option=o, dhw_scenario=d, metric=ADRIA.At([:fd_worst, :fd_best, :fd_net])] .= delta_tail_ratio(
                opt_fd, cf_fd; tail_fraction=tail_fraction
            )
        end
    end

    return aggreg, per_reef
end

"""
    plot_robustness(aggreg, N_seed_total, min_iv, rcp, horizon)

Aggregate 9-column line figure from `cvar_aggregation`'s `aggreg` YAXArray.
"""
function plot_robustness(aggreg, N_seed_total, min_iv, rcp, horizon)
    nseed_label = "1e$(round(Int, log10(N_seed_total)))"

    fig = Figure(size=(1300, 480))
    ax = Axis(fig[1, 1];
        xticks=(1:9, aggreg_metric_labels),
        xticklabelsize=10,
        ylabel="Relative performance against \ncounterfactual (no interv.)",
        title="Relative performance Δ/counterfactual - RCP $(rcp), N_seed $(nseed_label), min_locs $(min_iv), horizon $(horizon) yrs"
    )

    hlines!(ax, [0]; color=:black, linewidth=2)
    vlines!(ax, [3.5, 6.5]; color=:gray70, linewidth=1, linestyle=:dash)

    for (o_i, _) in enumerate(scenario_names)
        color = option_colors[o_i]
        for (d_i, _) in enumerate(dhw_scenarios)
            y = aggreg.data[o_i, d_i, :]
            lines!(ax, 1:9, y; color, linestyle=dhw_linestyles[d_i], linewidth=1.5)
            scatter!(ax, 1:9, y; color, marker=dhw_markers[d_i], markersize=10)
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
        joinpath(pd_config["plot_output_path"], "robustness_rcp$(rcp)_n$(nseed_label)_m$(min_iv)_h$(horizon).png"),
        fig; px_per_unit=2
    )
    @info "Saved robustness_rcp$(rcp)_n$(nseed_label)_m$(min_iv)_h$(horizon).png"
end

"""
    plot_boxplot(per_reef, N_seed_total, min_iv, rcp, horizon)

Per-reef boxplot (3 panels) from `cvar_aggregation`'s `per_reef` YAXArray, pooling all
reefs across DHW scenarios.
"""
function plot_boxplot(per_reef, N_seed_total, min_iv, rcp, horizon)
    nseed_label = "1e$(round(Int, log10(N_seed_total)))"

    fig = Figure(size=(1200, 450))
    for (m_i, metric) in enumerate(boxplot_metric_labels)
        ax = Axis(fig[1, m_i];
            xticks=(1:n_options, option_labels),
            xticklabelrotation=π / 4,
            ylabel=m_i == 1 ? "Relative difference vs counterfactual" : "",
            title="$(metric)"
        )
        hlines!(ax, [0]; color=:black, linewidth=2)

        for (o_i, _) in enumerate(scenario_names)
            # pool all reefs × all dhw_scenarios for this option
            y = vec(per_reef[option=o_i, metric=m_i])
            boxplot!(ax, fill(o_i, length(y)), y; color=option_colors[o_i])
        end
    end

    save(
        joinpath(pd_config["plot_output_path"], "robustness_boxplot_rcp$(rcp)_n$(nseed_label)_m$(min_iv)_h$(horizon).png"),
        fig; px_per_unit=2
    )
    @info "Saved robustness_boxplot_rcp$(rcp)_n$(nseed_label)_m$(min_iv)_h$(horizon).png"
end

# ── Analysis over different parameter sets  ──────────────────────────────────
# P10 / median / P90 of each net metric over all options × dhw_scenario, one point per
# parameter set (N_seed, RCP, min_iv_locations), whiskers = P10–P90.

summary_configs = [
    (N_seed=1e6, rcp="45", n_loc=200, horizon = 20),
    (N_seed=1e6, rcp="70", n_loc=200, horizon = 20),
    (N_seed=1e8, rcp="45", n_loc=200, horizon = 20),
    (N_seed=1e8, rcp="70", n_loc=200, horizon = 20),
    (N_seed=1e6, rcp="45", n_loc=100, horizon = 20),
    (N_seed=1e8, rcp="45", n_loc=100, horizon = 20)
]

net_syms = [:nyrs_net, :tac_net, :fd_net]
net_labels = ["Years >20% (net)", "Cum. cover (net)", "Cum. evenness (net)"]
net_colors = Makie.wong_colors()[1:length(net_syms)]

n_cfg = length(summary_configs)
config_labels = ["1e$(round(Int, log10(cfg.N_seed)))\nRCP$(cfg.rcp)\nn$(cfg.n_loc)" for cfg in summary_configs]
stat_syms = [:p10, :median, :p90]
# summary_stats: (config, net metric, stat) — P10 / median / P90 over options × DHW.
summary_stats = ADRIA.DataCube(
    zeros(n_cfg, length(net_syms), length(stat_syms));
    config=config_labels, metric=net_syms, stat=stat_syms
)

for (c_i, cfg) in enumerate(summary_configs)
    aggreg, per_reef = cvar_aggregation(
        rs, m_tac, fd_arr, loc_hab_area_km2,
        cfg.N_seed, cfg.n_loc, cfg.rcp, cfg.horizon; tail_fraction=tail_fraction
    )

    plot_robustness(aggreg, cfg.N_seed, cfg.n_loc, cfg.rcp, cfg.horizon)
    plot_boxplot(per_reef, cfg.N_seed, cfg.n_loc, cfg.rcp, cfg.horizon)

    for (m_i, m) in enumerate(net_syms)
        vals = vec(aggreg[metric=ADRIA.At(m)])  # over all option × dhw_scenario
        summary_stats[config=c_i, metric=ADRIA.At(m), stat=ADRIA.At(:p10)]    = quantile(vals, 0.10)
        summary_stats[config=c_i, metric=ADRIA.At(m), stat=ADRIA.At(:median)] = median(vals)
        summary_stats[config=c_i, metric=ADRIA.At(m), stat=ADRIA.At(:p90)]    = quantile(vals, 0.90)
    end
end

fig = Figure(size=(1100, 500))
ax = Axis(fig[1, 1];
    xticks=(1:n_cfg, config_labels),
    xticklabelsize=10,
    ylabel="Net relative performance vs counterfactual",
    title="Net metrics (P10 / median / P90 over options × DHW) — horizon 20 yrs"
)
hlines!(ax, [0]; color=:black, linewidth=1)

dodge = length(net_syms) == 1 ? [0.0] : collect(range(-0.2, 0.2; length=length(net_syms)))
for (m_i, m) in enumerate(net_syms)
    x = (1:n_cfg) .+ dodge[m_i]
    med = summary_stats[metric=ADRIA.At(m), stat=ADRIA.At(:median)].data
    lo = med .- summary_stats[metric=ADRIA.At(m), stat=ADRIA.At(:p10)].data
    hi = summary_stats[metric=ADRIA.At(m), stat=ADRIA.At(:p90)].data .- med
    errorbars!(ax, x, med, lo, hi; color=net_colors[m_i], whiskerwidth=8)
    scatter!(ax, x, med; color=net_colors[m_i], markersize=10)
end

Legend(fig[1, 2],
    [MarkerElement(; marker=:circle, color=net_colors[i]) for i in 1:length(net_syms)],
    net_labels; title="Net metric", framevisible=false
)

save(
    joinpath(pd_config["plot_output_path"], "robustness_net_summary_h$(horizon).png"),
    fig; px_per_unit=2
)
@info "Saved robustness_net_summary_h$(horizon).png"
