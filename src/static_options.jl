#=
Static seeding options: robustness + detail analysis.

One static option is applied for the whole seeding window. Scenarios are swept over a parameter
grid (N_seed budget, min_iv_locations, RCP, DHW); for each parameter set every option is compared
against a no-intervention counterfactual of the same DHW member.

Two analyses run on the single scenario table:
1. Robustness — CVaR aggregation of every option vs its counterfactual over a set of parameter
   sets (`summary_configs`), producing the aggregate / boxplot / net-summary figures.
2. Detail plots — for a chosen set of parameter sets (`static_plot_configs`), the per-option GBR
   maps (total seeds, Δ years-above-target, Δ cumulative cover), seeding-frequency histograms,
   per-option time-series, and the per-scenario seeding GIF.

Scenario generation follows the parameter-sweep layout (row = one option under one parameter set);
each row carries an `:option` identifier (-1 = counterfactual, 0 = unguided, 1:5 = guided option).
Analysis lives in `data_processing/{robustness,static_options}.jl`; plotting in
`visualization/{robustness,static_options}.jl`. This script only orchestrates.
=#

include("src/common.jl")
using GeoMakie, GraphMakie, CairoMakie, NaturalEarth
include("src/data_processing/robustness.jl")
include("src/visualization/robustness.jl")
include("src/data_processing/static_options.jl")
include("src/visualization/static_options.jl")

# ── Constants ─────────────────────────────────────────────────────────────────

rcps = ["45"]
dhw_scenarios = [2, 7, 10]
seed_years = 20
N_seed_totals = [1e6, 1e7, 1e8]
min_iv_locations_list = [100, 300, 500]
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
    timeframe=(2022, 2022 + seed_years + 2 + 5)
)
fix_common_parameters!(dom)

ADRIA.fix_factor!(dom;
    # Baseline only; every scenario row overrides min_iv_locations explicitly below.
    min_iv_locations=min_iv_locations_list[1],
    seed_years=seed_years,
    seeding_devices_per_m2=5,
    a_adapt=5.0,
    # Baseline only; every scenario row overrides N_seed_* explicitly below.
    N_seed_TA=N_seed_totals[1] * N_seed_weights.N_seed_TA,
    N_seed_CA=N_seed_totals[1] * N_seed_weights.N_seed_CA,
    N_seed_CNA=N_seed_totals[1] * N_seed_weights.N_seed_CNA,
    N_seed_SM=N_seed_totals[1] * N_seed_weights.N_seed_SM,
    N_seed_LM=N_seed_totals[1] * N_seed_weights.N_seed_LM
)

options = ADRIA.analysis.option_seed_preference(; include_weights=true)

# Assigned by name: both `options` and the model spec derive their criteria from
# `fieldnames(ADRIA.SeedCriteriaWeights)`, so this stays correct if criteria are added or
# removed (positional indexing into `options` silently breaks when they are).
criteria = ADRIA.component_params(ADRIA.model_spec(dom), "SeedCriteriaWeights").fieldname

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
        scens[row, criteria] = collect(option[criteria])
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

# ── Shared quantities ─────────────────────────────────────────────────────────

all_centroids = ADRIA.centroids(dom.loc_data)
dhw_model_names = dom.dhw_scens.properties["model_names"][dhw_scenarios]

# ── Robustness analysis (CVaR over parameter sets) ────────────────────────────

tail_number = 150

# Base quantities computed once (horizon-independent), consumed by `cvar_aggregation`.
n_locs = size(rs.seed_log, :locations)
m_tac, fd_arr, loc_hab_area_km2, seed_start = load_base_metrics(rs)

# Visualisation styling (shared across horizons)
option_colors = Makie.wong_colors()[1:n_options]
dhw_linestyles = [:solid, :dash, :dot]
dhw_markers = [:circle, :rect, :utriangle]
option_labels = string.(scenario_names)

# P10 / median / P90 of each net metric over all options × dhw_scenario, one point per
# parameter set (N_seed, RCP, min_iv_locations), whiskers = P10–P90.
summary_configs = [
    (N_seed=1e6, rcp="45", n_loc=100, horizon=20),
    (N_seed=1e6, rcp="45", n_loc=500, horizon=20),
    (N_seed=1e7, rcp="45", n_loc=500, horizon=20),
    (N_seed=1e7, rcp="45", n_loc=500, horizon=20),
    (N_seed=1e8, rcp="45", n_loc=100, horizon=20),
    (N_seed=1e8, rcp="45", n_loc=500, horizon=20)
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
        cfg.N_seed, cfg.n_loc, cfg.rcp, cfg.horizon; tail_number=tail_number
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

plot_net_summary(
    summary_stats, config_labels, net_syms, net_labels, net_colors, summary_configs[1].horizon
)

# ── Static-option detail plots (per selected parameter set) ───────────────────
# Each config is one parameter set plotted across all DHW members. `gif_dhw` selects which DHW
# member (index into `dhw_scenarios`) the per-scenario seeding GIF animates.

static_plot_configs = [
    (RCP="45", N_seed=1e6, min_iv=300, gif_dhw=1),
    (RCP="45", N_seed=1e7, min_iv=300, gif_dhw=1),
    (RCP="45", N_seed=1e8, min_iv=300, gif_dhw=1)
]

for cfg in static_plot_configs
    sel = select_static_scenarios(
        rs, cfg.N_seed, cfg.min_iv, cfg.RCP;
        N_seed_weights=N_seed_weights, dhw_scenarios=dhw_scenarios, scenario_names=scenario_names
    )
    seed = seeding_derived(rs, sel)
    perf = static_performance_metrics(rs, sel)
    tsm = scenario_timeseries_metrics(rs, seed.selected_locations)

    # Validation: seeding coverage across all intervention scenarios of this parameter set
    @info("Parameter set RCP $(cfg.RCP), N_seed $(cfg.N_seed), min_iv $(cfg.min_iv)")
    @info("  Total reefs  : $(size(dom.loc_data, 1))")
    @info("  Never seeded : $(sum(seed.never_seeded)) ($(round(100*mean(seed.never_seeded); digits=1))%)")
    @info("  Always seeded: $(sum(seed.always_seeded)) ($(round(100*mean(seed.always_seeded); digits=1))%)")

    plot_seeding_frequency_per_option(seed.seed_per_reef_per_ts_scen, sel, cfg)
    plot_seeding_frequency_all(seed.seed_per_reef_per_ts_scen, cfg)
    plot_total_seeds_map(seed.total_seeds, sel, cfg; seeds_dhw=1, colorrange=(0.0, 1e6))
    plot_nyrs_above_maps(perf.n_yrs_diff, seed.total_seeds, seed.active_mask, sel, cfg)
    plot_cum_tac_diff_maps(perf.cum_tac_diff, seed.total_seeds, seed.active_mask, sel, cfg)
    plot_metric_histograms(perf, sel, cfg)
    plot_options_timeseries(rs, tsm, sel, options, cfg)
    #animate_seeding_maps(rs, dom, seed.seed_per_reef_per_ts_scen, sel, cfg; gif_dhw=cfg.gif_dhw)
end
