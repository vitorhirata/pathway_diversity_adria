#=
Robustness analysis for pathway seeding.

Unlike `static_options.jl` (one static option applied for the whole
seeding window), here each scenario is a *pathway*: seeding runs for `seed_years` and the
option changes every `pd_frequency` years. Pathways are grouped by their **starting option**
(the option active at `seed_year_start`), and we show the range of per-reef performance a
starting option can lead to across all its downstream pathways.

Same three metrics as the static analysis, each compared against a no-intervention
counterfactual: years >20% coral cover, cumulative absolute cover, cumulative coral evenness.

Analysis lives in `robustness_data_processing.jl`; plotting in `robustness_visualisation.jl`.
This script only orchestrates: load results → compute metrics → make plots.
=#

include("src/common.jl")
using GeoMakie, GraphMakie, CairoMakie
include("src/data_processing/robustness.jl")
include("src/visualization/robustness.jl")

# ── Constants ─────────────────────────────────────────────

# These must match the pd_main.jl run that produced `intervention_path` below.
seed_years = 20
pd_frequency::Int64 = 5
seed_year_start = 2
N_seed_weights = (
    N_seed_TA=0.15, N_seed_CA=0.5, N_seed_CNA=0.0, N_seed_SM=0.35, N_seed_LM=0.0
)

# Measurement horizon (years from seed-start). Matches the static analysis' default horizon so
# the two scripts integrate over the same window and stay directly comparable.
horizon = seed_years

# ── Load intervention ResultSet (from a prior pd_main.jl run) ─────────────────

intervention_path = ""  # ADRIA results store produced by pd_main.jl
rs = ADRIA.load_results(intervention_path)

# ── Representative analysis parameters (the "main parameters") ────────────────

option_names = ADRIA.analysis.option_seed_preference().option_name
n_options = length(option_names)

# Parameters to analyse
tail_number = 150

# Parameter sets to plot: one set of figures A/B/C per entry. Each is a NamedTuple
# `(dhw, N_seed, n_locations)`. Mirrors pd_main.jl's `boxplot_configs`/`lockin_configs` style.
param_sets = [
    (dhw=7, N_seed=1e6, n_locations=100),
    (dhw=7, N_seed=1e7, n_locations=100),
    (dhw=7, N_seed=1e7, n_locations=500),
    (dhw=7, N_seed=1e8, n_locations=100),
    (dhw=7, N_seed=1e8, n_locations=500)
]

# ── Compute per-reef metrics (parameter-independent) ──────────────────────────

nyrs, ctac, cfd = reef_metrics(rs)
opt_metric_mats = (nyrs, ctac, cfd)

# Fixed colour per starting option, shared across all plots (first 5 Wong colors).
option_color_map = Dict(
    "heat_stress" => Makie.wong_colors()[1],
    "geographic_spread" => Makie.wong_colors()[2],
    "connectivity" => Makie.wong_colors()[3],
    "functional_diversity" => Makie.wong_colors()[4],
    "balanced" => Makie.wong_colors()[5]
)
option_colors = [option_color_map[string(o)] for o in option_names]
# Pretty option labels: "heat_stress" → "Heat stress"
option_labels = [uppercasefirst(replace(string(o), "_" => " ")) for o in option_names]

# ── Figures A/B/C — one set of figures per parameter set ──────────────────────
#
#   For each parameter set (dhw, N_seed, n_locations): group pathways by starting option, then
#   draw Figure A (per-reef delta boxplots), Figure B (CVaR tail-ratio ranges) and Figure C
#   (probability-weighted tail stats). Filenames encode the used dhw / N_seed / n_locations.

# Pathway adoption probabilities (Figure C), keyed by option_ts. From pd_main.jl's
# scenario_probabilities.csv.
prob_csv_path = ""

# Accumulate the plotted Figure B / Figure C data across parameter sets; written as one tidy CSV
# each after the loop. `_add_params` prepends the parameter-set id columns.
cvar_tables = DataFrame[]
weighted_tables = DataFrame[]
_add_params(df, ps) = insertcols(
    df, 1, :dhw => ps.dhw, :N_seed => round(Int, ps.N_seed), :n_locations => ps.n_locations
)

for ps in param_sets
    option_pathways = group_pathways_by_starting_option(
        rs, option_names;
        sel_dhw=ps.dhw, sel_N_seed=ps.N_seed, sel_n_locations=ps.n_locations,
        seed_year_start, seed_years, pd_frequency, N_seed_weights
    )
    cf_idx = find_counterfactual(rs, ps.dhw)
    cf_metric_full = (nyrs[:, cf_idx], ctac[:, cf_idx], cfd[:, cf_idx])

    # Figure A — per-reef delta boxplots pooled over pathways
    plot_pathways_boxplot(
        option_names, option_pathways, opt_metric_mats, cf_metric_full,
        option_colors, option_labels, ps
    )

    # Figure B — per starting option, median with P10–P90 whiskers of per-pathway ratios,
    # across the 9 metric variants (GBR/worst/best × 3 metrics).
    cvar_df = pathways_cvar_ranges(
        option_names, option_pathways, opt_metric_mats, cf_metric_full;
        tail_number, seed_years
    )
    plot_pathways_cvar(cvar_df, option_names, option_colors, option_labels, ps)
    push!(cvar_tables, _add_params(cvar_df, ps))

    # Figure C — probability-weighted tail stats. Weight the per-pathway scalars by pathway
    # adoption probability and report P10 / Median / Mean / P90, GBR/worst/best kept separate.
    prob_map = load_pathway_probabilities(
        prob_csv_path; sel_N_seed=ps.N_seed, sel_dhw=ps.dhw, sel_n_locations=ps.n_locations
    )
    # (metric name, intervention matrix, counterfactual column, reference constant)
    weighted_metrics = [
        ("Years >20% coral cover", nyrs, cf_metric_full[1], float(seed_years)),
        ("Cumulative cover", ctac, cf_metric_full[2], mean(cf_metric_full[2])),
        ("Cumulative evenness", cfd, cf_metric_full[3], mean(cf_metric_full[3]))
    ]
    weighted_tail_stats = compute_weighted_tail_stats(
        option_names, option_pathways, weighted_metrics, prob_map, rs; tail_number=tail_number
    )
    plot_pathways_weighted(
        weighted_tail_stats, option_names, option_colors, option_labels, ps;
        point_stat=:median
    )
    # Tidy Figure C rows: machine metric label from metric_idx, keep median/mean/P10/P90.
    # `select` is qualified: DataFrames and another loaded package (YAXArrays) both export it,
    weighted_out = DataFrames.select(
        weighted_tail_stats,
        :start_option, :metric_idx,
        :metric_idx => DataFrames.ByRow(i -> cvar_metric_labels[i]) => :metric,
        :median, :mean, :p10, :p90
    )
    push!(weighted_tables, _add_params(weighted_out, ps))
end

# ── Save plotted Figure B / Figure C data (tidy, one row per option × metric variant) ─────────
CSV.write(
    joinpath(pd_config["plot_output_path"], "robustness_scores.csv"),
    vcat(cvar_tables...)
)
CSV.write(
    joinpath(pd_config["plot_output_path"], "robustness_scores_weighted_prob.csv"),
    vcat(weighted_tables...)
)
@info "Saved robustness_scores.csv and robustness_scores_weighted_prob.csv"

# ── Aggregated robustness across DHW scenarios and parameter sets ─────────────
#
#   Per pathway, collapse the 6 tail ratios into one scalar (mean of the 6 bottom/top ratios) and
#   keep its worst (smallest) value across DHW scenarios. Then, for each (starting option ×
#   parameter set), report P10 / median / P90 of that worst-case robustness over pathways.

dhw_scenarios = [2, 7, 10]  # all DHW members present in the run
param_sets_nl = detect_param_sets(rs; N_seed_weights)  # (N_seed, n_locations), worst over DHW

rob_df = worst_dhw_robustness(
    rs, option_names, opt_metric_mats;
    dhw_scenarios, param_sets=param_sets_nl,
    seed_year_start, seed_years, pd_frequency, N_seed_weights, tail_number
)
CSV.write(joinpath(pd_config["plot_output_path"], "robustness.csv"), rob_df)
# Combine two simulations by stacking rows
#rob_df = vcat(rob_df, rob_df2)

plot_robustness_param_scatter(rob_df, option_names, option_colors, option_labels)

# ── Top robustness pathways per starting option and parameter set ─────────────
#
#   For each (starting option × parameter set), rank the per-block option sequences by their
#   worst-over-DHW robustness (min over the DHW scenarios each sequence occurs in) and keep the
#   top 10. Sequences are the compact 4-block pathways (e.g. "heat_stress > balanced > ...").
#   One CSV per (N_seed, n_locations) parameter set.

for ps in param_sets_nl
    top_pathways = top_robustness_pathways(
        rs, option_names, opt_metric_mats;
        dhw_scenarios, param_sets=[ps],
        seed_year_start, seed_years, pd_frequency, N_seed_weights, tail_number
    )
    fname = "top_robustness_pathways_n$(_sci(ps.N_seed))_l$(ps.n_locations).csv"
    CSV.write(joinpath(pd_config["plot_output_path"], fname), top_pathways)
    @info "Saved $(fname)"
end

# ── Figure D — robustness vs pathway diversity, facetted by parameter set ─────
#
#   Load a pathway-diversity table (same format as pd_main.jl's `options`:
#   option_name, pathway_diversity, N_seed, dhw_scenario, n_locations, rcp) and, per parameter
#   set (panel), scatter worst-DHW median robustness (y) against worst-DHW pathway diversity (x),
#   coloured by starting option. Panels are gridded: columns = n_locations, rows = N_seed.

pathway_diversity_csv_path = ""  # pathway_diversity.csv written by pd_main.jl

pd_df = CSV.read(pathway_diversity_csv_path, DataFrame)
rob_div_df = join_robustness_diversity(rob_df, pd_df)

# Rank agreement between pathway diversity and robustness across starting options. The stratified
# (blocked) Kendall τ-b compares options only within a shared N_seed × n_locations parameter set,
# summing pair counts across panels — holding the parameter-set gradient fixed (descriptive; n is
# small per panel, so no significance test).
tau_agreement = kendall_tau_diversity_robustness(rob_div_df)
@info "Stratified (blocked) Kendall τ-b, robustness vs pathway diversity" τ = tau_agreement.stratified

plot_robustness_vs_diversity(rob_div_df, option_names, option_colors, option_labels)
