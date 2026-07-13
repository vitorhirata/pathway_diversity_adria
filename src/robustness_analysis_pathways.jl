#=
Robustness analysis for pathway seeding.

Unlike `robustness_analysis_static_options.jl` (one static option applied for the whole
seeding window), here each scenario is a *pathway*: seeding runs for `seed_years` and the
option changes every `pd_frequency` years. Pathways are grouped by their **starting option**
(the option active at `seed_year_start`), and we show the range of per-reef performance a
starting option can lead to across all its downstream pathways.

Same three metrics as the static analysis, each compared against a no-intervention
counterfactual: years >20% coral cover, cumulative absolute cover, cumulative coral evenness.
=#

include("src/common.jl")
using GeoMakie, GraphMakie, CairoMakie

# ── Constants ─────────────────────────────────────────────

rcps = ["26", "45", "70"]
dhw_scenarios = [5, 7, 11]
seed_years = 20
pd_frequency::Int64 = 5
seed_year_start = 2
N_seed_weights = (
    N_seed_TA=0.15, N_seed_CA=0.5, N_seed_CNA=0.0, N_seed_SM=0.35, N_seed_LM=0.0
)

# ── Domain setup (same fixed factors and timeframe as pd_main.jl) ─────────────

dom = ADRIA.load_domain(
    pd_config["domain_path"], rcps[1];
    calib_params_fn=pd_config["coral_param_path"],
    timeframe=(2022, 2022 + seed_years + 2 + 5)
)
ms = ADRIA.model_spec(dom)

ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "FogCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "MCCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "Coral").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "GrowthAcceleration").fieldname)

ADRIA.fix_factor!(dom;
    seed_year_start=seed_year_start,
    seed_years=seed_years,
    seed_deployment_freq=1,
    seeding_devices_per_m2=5,
    seed_strategy=1,
    a_adapt=5.0,
    a_adapt_ref=5,
    plan_horizon=5.0,
    guided=1,
    depth_min=2.0,
    depth_offset=25.0,
    fogging=0.0,
    SRM=0.0,
    N_mc_settlers=0,
    wave_scenario=1
)

# ── Build + run counterfactual ResultSet (one per dhw × rcp, N_seed = 0) ──────

scens_cf = repeat(ADRIA.sample(dom, 2)[1:1, :], length(dhw_scenarios))
for (i, dhw) in enumerate(dhw_scenarios)
    scens_cf[i, :dhw_scenario] = dhw
    scens_cf[i, :guided] = 0
    scens_cf[i, :N_seed_TA] = 0
    scens_cf[i, :N_seed_CA] = 0
    scens_cf[i, :N_seed_CNA] = 0
    scens_cf[i, :N_seed_SM] = 0
    scens_cf[i, :N_seed_LM] = 0
end

rs_cf = ADRIA.run_scenarios(dom, scens_cf, rcps)

# rs_cf = ADRIA.load_results()

# ── Load intervention ResultSet (from a prior pd_main.jl run) ─────────────────

intervention_path = ""  # ADRIA results store produced by pd_main.jl
rs = ADRIA.load_results(intervention_path)

# ── Representative analysis parameters (the "main parameters") ────────────────

sel_rcp = 45
sel_dhw = 5
sel_N_seed = 1e7
sel_n_locations = 200
tail_fraction = 0.03

option_names = ADRIA.analysis.option_seed_preference().option_name
n_options = length(option_names)

metric_labels = [
    "Years >20% coral cover\n(worst reefs)", "Years >20% coral cover\n(best reefs)",
    "Cumulative cover\n(worst reefs)", "Cumulative cover\n(best reefs)",
    "Cumulative evenness\n(worst reefs)", "Cumulative evenness\n(best reefs)"
]

"""
Select bottom and top `tail_fraction` of reefs by delta (opt − cf), then compute
(mean(opt[group]) − mean(cf[group])) / mean(cf[group]) for each group.
Returns (bottom_ratio, top_ratio).
"""
function delta_tail_ratio(
    opt::AbstractVector, cf::AbstractVector; tail_fraction::Float64=0.05, norm::Real = 0.0
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

# ── Per-reef metrics (same code path for intervention and counterfactual) ─────

"""
Compute the three per-reef metrics for a ResultSet.
Returns (n_yrs_above, cum_tac, cum_fd), each a (n_locs, n_scenarios) matrix.
"""
function reef_metrics(result_set)
    n_locs = size(result_set.seed_log, :locations)
    n_scens = nrow(result_set.inputs)

    m_tac = Array(ADRIA.metrics.total_absolute_cover(result_set)) .* 1e-6
    loc_hab_area_km2 = result_set.loc_area .* result_set.loc_max_coral_cover .* 1e-6
    fd_data = ADRIA.metrics.coral_evenness(result_set)

    n_yrs_above = Array{Float64}(undef, n_locs, n_scens)
    for s in 1:n_scens, l in 1:n_locs
        thr = 0.20 * loc_hab_area_km2[l]
        n_yrs_above[l, s] = Float64(count(m_tac[:, l, s] .>= thr))
    end

    cum_tac = dropdims(sum(m_tac; dims=1); dims=1)
    cum_fd = Array(dropdims(sum(fd_data; dims=:timesteps); dims=:timesteps))

    return n_yrs_above, cum_tac, cum_fd
end

nyrs, ctac, cfd = reef_metrics(rs)
nyrs_cf, ctac_cf, cfd_cf = reef_metrics(rs_cf)

# ── Group pathways by starting option ─────────────────────────────────────────

max_time = size(rs.seed_log, :timesteps)
starting_option(option_ts) = ADRIA.analysis.decode_option_ts(
    option_ts, seed_year_start, seed_years, pd_frequency, max_time
)[seed_year_start]

# Intervention scenarios matching the representative combo
combo_mask =
    rs.inputs.N_seed_CA .== sel_N_seed * N_seed_weights.N_seed_CA .&&
    rs.inputs.N_seed_TA .== sel_N_seed * N_seed_weights.N_seed_TA .&&
    rs.inputs.min_iv_locations .== sel_n_locations .&&
    rs.inputs.dhw_scenario .== sel_dhw .&&
    rs.inputs.RCP .== sel_rcp
combo_idxs = findall(combo_mask)

# Scenario indices grouped by starting option
option_pathways = Dict(
    option => [i for i in combo_idxs if starting_option(rs.inputs.option_ts[i]) == option]
    for option in option_names
)

# Matching counterfactual column (same dhw + rcp)
cf_idx = findfirst(
    (rs_cf.inputs.dhw_scenario .== sel_dhw) .& (rs_cf.inputs.RCP .== sel_rcp)
)

# Per-reef metric matrices (intervention) and matching counterfactual columns
opt_metric_mats = (nyrs, ctac, cfd)
cf_metric_full = (nyrs_cf[:, cf_idx], ctac_cf[:, cf_idx], cfd_cf[:, cf_idx])

option_colors = Makie.wong_colors()[1:n_options]
option_labels = string.(option_names)

# ── Figure A — boxplot per starting option (per-reef deltas pooled over paths) ─

metric_titles = [
    "Years >20% coral cover", "Cumulative cover", "Cumulative evenness"
]

fig = Figure(size=(1200, 450))
for (m_i, title) in enumerate(metric_titles)
    ax = Axis(fig[1, m_i];
        xticks=(1:n_options, option_labels),
        xticklabelrotation=π / 4,
        ylabel=m_i == 1 ? "Δ vs counterfactual (per reef)" : "",
        title=title
    )
    hlines!(ax, [0]; color=:black, linewidth=2)

    for (o_i, option) in enumerate(option_names)
        idxs = option_pathways[option]
        isempty(idxs) && continue
        # Pool per-reef deltas over all reefs × all pathways in this group
        y = Float64[]
        for s in idxs
            append!(y, opt_metric_mats[m_i][:, s] .- cf_metric_full[m_i])
        end
        boxplot!(ax, fill(o_i, length(y)), y; color=option_colors[o_i])
    end
end
save(
    joinpath(
        pd_config["plot_output_path"],
        "robustness_pathways_boxplot_rcp$(sel_rcp)_dhw$(sel_dhw).png"
    ),
    fig; px_per_unit=2
)
@info "Saved robustness_pathways_boxplot_rcp$(sel_rcp)_dhw$(sel_dhw).png"

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
        b1, t1 = delta_tail_ratio(opt_metric_mats[1][:, s], cf_metric_full[1]; tail_fraction=tail_fraction, norm=seed_years)
        b2, t2 = delta_tail_ratio(opt_metric_mats[2][:, s], cf_metric_full[2]; tail_fraction=tail_fraction)
        b3, t3 = delta_tail_ratio(opt_metric_mats[3][:, s], cf_metric_full[3]; tail_fraction=tail_fraction)
        ratios[p, :] = [b1, t1, b2, t2, b3, t3]
    end

    for m in 1:6
        med[o_i, m] = median(ratios[:, m])
        lo[o_i, m] = minimum(ratios[:, m])
        hi[o_i, m] = maximum(ratios[:, m])
    end
end

fig = Figure(size=(1100, 480))
ax = Axis(fig[1, 1];
    xticks=(1:6, metric_labels),
    xticklabelsize=11,
    ylabel="Relative performance against \ncounterfactual (no interv.)",
    title="Pathway robustness by starting option — RCP $(sel_rcp), dhw $(sel_dhw)"
)
hlines!(ax, [0]; color=:black, linewidth=2)
vlines!(ax, [2.5, 4.5]; color=:gray70, linewidth=1, linestyle=:dash)

n_dodge = n_options
dodge_width = 0.6
for (o_i, _) in enumerate(option_names)
    all(isnan, med[o_i, :]) && continue
    color = option_colors[o_i]
    offset = (o_i - (n_dodge + 1) / 2) * (dodge_width / n_dodge)
    x = (1:6) .+ offset
    scatter!(ax, x, med[o_i, :]; color, markersize=10)
    errorbars!(
        ax, x, med[o_i, :], med[o_i, :] .- lo[o_i, :], hi[o_i, :] .- med[o_i, :];
        color, whiskerwidth=6
    )
end

Legend(fig[1, 2],
    [MarkerElement(; marker=:circle, color=option_colors[i]) for i in 1:n_options],
    option_labels;
    title="Starting option", framevisible=false
)

save(
    joinpath(
        pd_config["plot_output_path"],
        "robustness_pathways_rcp$(sel_rcp)_dhw$(sel_dhw).png"
    ),
    fig; px_per_unit=2
)
@info "Saved robustness_pathways_rcp$(sel_rcp)_dhw$(sel_dhw).png"

# ── Probability-weighted tail statistics per starting option ──────────────────
#
#   Weight the per-pathway scalars by pathway adoption probability and
#     report P10 / Median / Mean / P90 of the weighted distribution, top and bottom kept
#     separate. Done per starting option.

"""
    weighted_quantile(values, weights, q)

Step / inverse-CDF weighted quantile (no interpolation): smallest value whose cumulative
weight (values sorted ascending) reaches `q`.
"""
function weighted_quantile(values::AbstractVector, weights::AbstractVector, q::Real)
    order = sortperm(values)
    vs = values[order]
    ws = weights[order]
    F = cumsum(ws)
    i = findfirst(f -> f >= q - 1e-9, F)
    return vs[i]
end

"""
    pathway_tail_stats(M, w) -> (p10, median, mean, p90)

Weighted stats of a per-pathway scalar vector `M` with weights `w`. `mean` uses the unsorted
vectors directly; the quantiles use `weighted_quantile`.
"""
function pathway_tail_stats(M::AbstractVector, w::AbstractVector)
    return (
        p10=weighted_quantile(M, w, 0.10),
        median=weighted_quantile(M, w, 0.50),
        mean=sum(w .* M),
        p90=weighted_quantile(M, w, 0.90)
    )
end

"""
    aggregate_metric(intervention, counterfactual, reference, probs; tail_frac=0.05)

`intervention` is `nreefs × npathways`, `counterfactual` is length `nreefs`, `probs` is length
`npathways`. Returns `(top, bottom, n_eff)` where `top`/`bottom` are `pathway_tail_stats`
named tuples for the best/worst reef tails and `n_eff` is the Kish effective sample size.
"""
function aggregate_metric(
    intervention::AbstractMatrix, counterfactual::AbstractVector, reference::Real,
    probs::AbstractVector; tail_frac::Float64=0.05
)
    nreefs, npaths = size(intervention)
    @assert nreefs == length(counterfactual)
    @assert npaths == length(probs)
    @assert all(probs .>= 0)
    @assert isapprox(sum(probs), 1.0; atol=1e-5)

    Mplus = Vector{Float64}(undef, npaths)   # top tail (best reefs)
    Mminus = Vector{Float64}(undef, npaths)  # bottom tail (worst reefs)
    for p in 1:npaths
        Mminus[p], Mplus[p] = delta_tail_ratio(
            intervention[:, p], counterfactual; tail_fraction=tail_frac, norm=reference
        )
    end

    return (
        top=pathway_tail_stats(Mplus, probs),
        bottom=pathway_tail_stats(Mminus, probs),
        n_eff=1 / sum(probs .^ 2)
    )
end

# Pathway adoption probabilities for the representative combo, keyed by option_ts
prob_csv_path = ""

prob_df = CSV.read(prob_csv_path, DataFrame)
prob_df = prob_df[
    (prob_df.N_seed .== round(Int, sel_N_seed)) .&
    (prob_df.dhw_scenario .== sel_dhw) .&
    (prob_df.n_locations .== sel_n_locations) .&
    (prob_df.rcp .== sel_rcp), :]
prob_map = Dict(row.option_ts => row.probability for row in eachrow(prob_df))

# (metric name, intervention matrix, counterfactual column, reference constant)
weighted_metrics = [
    ("Years >20% coral cover", nyrs, cf_metric_full[1], float(seed_years)),
    ("Cumulative cover", ctac, cf_metric_full[2], mean(cf_metric_full[2])),
    ("Cumulative evenness", cfd, cf_metric_full[3], mean(cf_metric_full[3]))
]

# Tidy per (starting option × metric × tail) table. metric_idx 1–6 matches `metric_labels`
# ordering (worst/best interleaved per metric), so this feeds a Figure-B-style plot directly.
weighted_tail_stats = DataFrame(;
    start_option=String[], metric=String[], metric_idx=Int[], tail=String[],
    p10=Float64[], median=Float64[], mean=Float64[], p90=Float64[], n_eff=Float64[]
)

for option in option_names
    idxs = option_pathways[option]
    isempty(idxs) && continue
    probs = [prob_map[rs.inputs.option_ts[s]] for s in idxs]

    for (m, (name, mat, cf, ref)) in enumerate(weighted_metrics)
        res = aggregate_metric(mat[:, idxs], cf, ref, probs; tail_frac=tail_fraction)
        for (tail_name, stats, metric_idx) in (
            ("worst", res.bottom, 2 * (m - 1) + 1), ("best", res.top, 2 * m)
        )
            push!(weighted_tail_stats, (
                string(option), name, metric_idx, tail_name,
                stats.p10, stats.median, stats.mean, stats.p90, res.n_eff
            ))
        end
    end
end

@info "Computed weighted_tail_stats" nrow(weighted_tail_stats)

# ── Figure C — probability-weighted tail stats per starting option ────────────
# Same layout as Figure B, but points come from the weighted distribution: point = median
# (or mean), whiskers span P10–P90. Set `point_stat` to :median or :mean.

point_stat = :median  # :median or :mean

fig = Figure(size=(1100, 480))
ax = Axis(fig[1, 1];
    xticks=(1:6, metric_labels),
    xticklabelsize=11,
    ylabel="Prob.-weighted performance against \ncounterfactual (no interv.)",
    title="Weighted pathway robustness ($(point_stat), P10–P90) — RCP $(sel_rcp), dhw $(sel_dhw)"
)
hlines!(ax, [0]; color=:black, linewidth=2)
vlines!(ax, [2.5, 4.5]; color=:gray70, linewidth=1, linestyle=:dash)

n_dodge = n_options
dodge_width = 0.6
for (o_i, option) in enumerate(option_names)
    df = weighted_tail_stats[weighted_tail_stats.start_option .== string(option), :]
    isempty(df) && continue
    sort!(df, :metric_idx)

    color = option_colors[o_i]
    offset = (o_i - (n_dodge + 1) / 2) * (dodge_width / n_dodge)
    x = df.metric_idx .+ offset
    pt = df[!, point_stat]
    scatter!(ax, x, pt; color, markersize=10)
    errorbars!(ax, x, pt, pt .- df.p10, df.p90 .- pt; color, whiskerwidth=6)
end

Legend(fig[1, 2],
    [MarkerElement(; marker=:circle, color=option_colors[i]) for i in 1:n_options],
    option_labels;
    title="Starting option", framevisible=false
)

save(
    joinpath(
        pd_config["plot_output_path"],
        "robustness_pathways_weighted_$(point_stat)_rcp$(sel_rcp)_dhw$(sel_dhw).png"
    ),
    fig; px_per_unit=2
)
@info "Saved robustness_pathways_weighted_$(point_stat)_rcp$(sel_rcp)_dhw$(sel_dhw).png"
