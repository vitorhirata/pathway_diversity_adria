#=
Robustness analysis for pathway seeding.

Unlike `robustness_analysis_static_options.jl` (one static option applied for the whole
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
include("src/helpers/robustness_data_processing.jl")
include("src/helpers/robustness_visualisation.jl")

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
tail_fraction = 0.03
sel_rcp = 45
sel_dhw = 7
sel_N_seed = 1e6
sel_n_locations = 200

# ── Compute per-reef metrics and group pathways ───────────────────────────────

nyrs, ctac, cfd = reef_metrics(rs)

option_pathways = group_pathways_by_starting_option(
    rs, option_names;
    sel_rcp, sel_dhw, sel_N_seed, sel_n_locations,
    seed_year_start, seed_years, pd_frequency, N_seed_weights
)
cf_idx = find_counterfactual(rs, sel_dhw, sel_rcp)

# Per-reef metric matrices (intervention) and matching counterfactual columns
opt_metric_mats = (nyrs, ctac, cfd)
cf_metric_full = (nyrs[:, cf_idx], ctac[:, cf_idx], cfd[:, cf_idx])

option_colors = Makie.wong_colors()[1:n_options]
option_labels = string.(option_names)

# ── Figure A — boxplot per starting option (per-reef deltas pooled over paths) ─

plot_pathways_boxplot(
    option_names, option_pathways, opt_metric_mats, cf_metric_full,
    option_colors, option_labels, sel_rcp, sel_dhw
)

# ── Figure B — CVaR range with error bars per starting option ─────────────────

# Per starting option: median with min/max whiskers of per-pathway CVaR ratios,
# across the 6 metric variants (worst/best × 3 metrics).
med = fill(NaN, n_options, 6)
lo = fill(NaN, n_options, 6)
hi = fill(NaN, n_options, 6)

for (o_i, option) in enumerate(option_names)
    idxs = option_pathways[option]
    isempty(idxs) && continue

    # per-pathway CVaR ratios: rows = pathways, cols = 6 metric variants
    ratios = Array{Float64}(undef, length(idxs), 6)
    for (p, s) in enumerate(idxs)
        b1, t1, _ = delta_tail_ratio(opt_metric_mats[1][:, s], cf_metric_full[1]; tail_fraction=tail_fraction, norm=seed_years)
        b2, t2, _ = delta_tail_ratio(opt_metric_mats[2][:, s], cf_metric_full[2]; tail_fraction=tail_fraction)
        b3, t3, _ = delta_tail_ratio(opt_metric_mats[3][:, s], cf_metric_full[3]; tail_fraction=tail_fraction)
        ratios[p, :] = [b1, t1, b2, t2, b3, t3]
    end

    for m in 1:6
        med[o_i, m] = median(ratios[:, m])
        lo[o_i, m] = minimum(ratios[:, m])
        hi[o_i, m] = maximum(ratios[:, m])
    end
end

plot_pathways_cvar(med, lo, hi, option_names, option_colors, option_labels, sel_rcp, sel_dhw)

# ── Aggregated robustness across DHW scenarios and parameter sets ─────────────
#
#   Per pathway, collapse the 6 tail ratios into one scalar (mean over metrics of top − |bottom|).
#   Then, for each (starting option × parameter set), take the median robustness over pathways per
#   DHW and keep the worst (smallest-median) DHW; report that DHW's median/min/max over pathways.

dhw_scenarios = [2, 7, 10]  # all DHW members present in the run
param_sets = detect_param_sets(rs, sel_rcp; N_seed_weights)

rob_df = worst_dhw_robustness(
    rs, option_names, opt_metric_mats;
    dhw_scenarios, param_sets, sel_rcp,
    seed_year_start, seed_years, pd_frequency, N_seed_weights, tail_fraction
)

plot_robustness_param_scatter(rob_df, option_names, option_colors, option_labels, sel_rcp)

# ── Top robustness pathways per starting option and parameter set ─────────────
#
#   For each (starting option × parameter set), rank the per-block option sequences by their
#   worst-over-DHW robustness (min over the DHW scenarios each sequence occurs in) and keep the
#   top 10. Sequences are the compact 4-block pathways (e.g. "heat_stress > balanced > ...").

top_pathways = top_robustness_pathways(
    rs, option_names, opt_metric_mats;
    dhw_scenarios, param_sets, sel_rcp,
    seed_year_start, seed_years, pd_frequency, N_seed_weights, tail_fraction
)
CSV.write(
    joinpath(pd_config["plot_output_path"], "top_robustness_pathways_rcp$(sel_rcp).csv"),
    top_pathways
)
@info "Saved top_robustness_pathways_rcp$(sel_rcp).csv"

# ── Probability-weighted tail statistics per starting option ──────────────────
#
#   Weight the per-pathway scalars by pathway adoption probability and
#     report P10 / Median / Mean / P90 of the weighted distribution, top and bottom kept
#     separate. Done per starting option.

# Pathway adoption probabilities for the representative combo, keyed by option_ts
prob_csv_path = ""

prob_map = load_pathway_probabilities(
    prob_csv_path; sel_N_seed, sel_dhw, sel_n_locations, sel_rcp
)

# (metric name, intervention matrix, counterfactual column, reference constant)
weighted_metrics = [
    ("Years >20% coral cover", nyrs, cf_metric_full[1], float(seed_years)),
    ("Cumulative cover", ctac, cf_metric_full[2], mean(cf_metric_full[2])),
    ("Cumulative evenness", cfd, cf_metric_full[3], mean(cf_metric_full[3]))
]

weighted_tail_stats = compute_weighted_tail_stats(
    option_names, option_pathways, weighted_metrics, prob_map, rs; tail_fraction=tail_fraction
)
@info "Computed weighted_tail_stats" nrow(weighted_tail_stats)

# ── Figure C — probability-weighted tail stats per starting option ────────────

plot_pathways_weighted(
    weighted_tail_stats, option_names, option_colors, option_labels, sel_rcp, sel_dhw;
    point_stat=:median
)

# ── Figure D — robustness vs pathway diversity, facetted by parameter set ─────
#
#   Load a pathway-diversity table (same format as pd_main.jl's `options`:
#   option_name, pathway_diversity, N_seed, dhw_scenario, n_locations, rcp) and, per parameter
#   set (panel), scatter worst-DHW median robustness (y) against worst-DHW pathway diversity (x),
#   coloured by starting option. Panels are gridded: columns = n_locations, rows = N_seed.

pathway_diversity_csv_path = ""  # pathway_diversity.csv written by pd_main.jl

pd_df = CSV.read(pathway_diversity_csv_path, DataFrame)
rob_div_df = join_robustness_diversity(rob_df, pd_df, sel_rcp)

plot_robustness_vs_diversity(rob_div_df, option_names, option_colors, option_labels, sel_rcp)
