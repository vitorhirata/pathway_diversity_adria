#=
Data processing for the robustness analyses.

Shared analysis / metric-computation code for both `robustness_analysis_pathways.jl` (pathways
grouped by starting option) and `robustness_analysis_static_options.jl` (one static option for
the whole seeding window). No plotting lives here — see `robustness_visualisation.jl`.

This file is *included by* the main scripts and assumes `include("src/common.jl")` has already
run (it uses `ADRIA`, `Statistics`' `mean`/`median`/`quantile`, `CSV`, `DataFrames`). Several
functions read analysis globals defined by the including script (`scenario_names`,
`dhw_scenarios`, `n_options`, `n_dhw`, `n_locs`, `N_seed_weights`, `seed_start`).
=#

# ── CVaR tail ratios ──────────────────────────────────────────────────────────

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
    opt::AbstractVector, cf::AbstractVector; tail_fraction::Float64=0.03, norm::Real=0.0
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

# ── Base metrics and windowing ────────────────────────────────────────────────

"""
    load_base_metrics(result_set) -> (m_tac, fd_arr, loc_hab_area_km2, seed_start)

Materialise the horizon-independent base quantities shared by both robustness analyses.
`readcubedata` keeps the metrics as YAXArrays (dims preserved) so named slicing stays cheap.
"""
function load_base_metrics(result_set)
    m_tac = ADRIA.readcubedata(ADRIA.metrics.total_absolute_cover(result_set)) .* 1e-6  # km²
    fd_arr = ADRIA.readcubedata(ADRIA.metrics.coral_evenness(result_set))
    loc_hab_area_km2 = result_set.loc_area .* result_set.loc_max_coral_cover .* 1e-6
    seed_start = Int(result_set.inputs.seed_year_start[1])
    return m_tac, fd_arr, loc_hab_area_km2, seed_start
end

"""Per-location count of timesteps in `window` where scenario `s`'s cover ≥ 20% habitable area."""
function window_nyrs(m_tac, loc_hab_area_km2, window, s)
    n_locs = size(m_tac, :locations)
    return Float64[
        count(m_tac[timesteps=window, locations=l, scenarios=s].data .>= 0.20 * loc_hab_area_km2[l])
        for l in 1:n_locs
    ]
end

"""Per-location cumulative absolute cover (km²·yr) over `window` for scenario `s`."""
window_tac(m_tac, window, s) = vec(sum(m_tac[timesteps=window, scenarios=s].data; dims=1))

"""Per-location cumulative coral evenness over `window` for scenario `s`."""
window_fd(fd_arr, window, s) = vec(sum(fd_arr[timesteps=window, scenarios=s].data; dims=1))

"""
    reef_metrics(result_set; horizon) -> (n_yrs_above, cum_tac, cum_fd)

Three per-reef metrics for a ResultSet over a `horizon`-year window anchored at seed-start, each
a `(n_locs, n_scenarios)` matrix. Used by the pathways analysis, which then compares columns
against a counterfactual column.
"""
function reef_metrics(result_set; horizon::Int=horizon)
    m_tac, fd_arr, loc_hab_area_km2, seed_start = load_base_metrics(result_set)
    n_locs = size(result_set.seed_log, :locations)
    n_scens = nrow(result_set.inputs)

    @assert seed_start + horizon - 1 <= size(m_tac, :timesteps) "Run too short for horizon $(horizon)."
    window = seed_start:(seed_start + horizon - 1)

    n_yrs_above = Array{Float64}(undef, n_locs, n_scens)
    cum_tac = Array{Float64}(undef, n_locs, n_scens)
    cum_fd = Array{Float64}(undef, n_locs, n_scens)
    for s in 1:n_scens
        n_yrs_above[:, s] = window_nyrs(m_tac, loc_hab_area_km2, window, s)
        cum_tac[:, s] = window_tac(m_tac, window, s)
        cum_fd[:, s] = window_fd(fd_arr, window, s)
    end

    return n_yrs_above, cum_tac, cum_fd
end

# ── Static-option CVaR aggregation ────────────────────────────────────────────

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

Reads analysis globals `scenario_names`, `dhw_scenarios`, `n_options`, `n_dhw`, `n_locs`,
`N_seed_weights`, `seed_start` from the including script.
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

        cf_nyrs = window_nyrs(m_tac, loc_hab_area_km2, window, cf_s)
        cf_tac  = window_tac(m_tac, window, cf_s)
        cf_fd   = window_fd(fd_arr, window, cf_s)

        for (o, _) in enumerate(scenario_names)
            # option id: 1:5 for guided options, 0 for the unguided (6th) column.
            opt_id = o < n_options ? o : 0
            s = findfirst(seeded_mask .& dhw_mask .& (rs.inputs.option .== opt_id))
            opt_nyrs = window_nyrs(m_tac, loc_hab_area_km2, window, s)
            opt_tac  = window_tac(m_tac, window, s)
            opt_fd   = window_fd(fd_arr, window, s)

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

# ── Pathway grouping (pathways analysis) ──────────────────────────────────────

"""
    compact_pathway(option_ts, seed_year_start, seed_years, pd_frequency, max_time) -> Vector

The per-block option sequence of a pathway: one option per `pd_frequency`-year decision block
within the seeding window, dropping the `:nothing` padding outside `seed_years`. Length is
`seed_years ÷ pd_frequency` (e.g. `[:heat_stress, :balanced, :functional_diversity, :heat_stress]`).
"""
function compact_pathway(option_ts, seed_year_start, seed_years, pd_frequency, max_time)
    number_changes = seed_years ÷ pd_frequency
    decision_steps = [seed_year_start + (k - 1) * pd_frequency for k in 1:number_changes]
    return ADRIA.analysis.decode_option_ts(
        option_ts, seed_year_start, seed_years, pd_frequency, max_time
    )[decision_steps]
end

"""
    group_pathways_by_starting_option(rs, option_names; sel_rcp, sel_dhw, sel_N_seed,
        sel_n_locations, seed_year_start, seed_years, pd_frequency, N_seed_weights)

Scenario indices of the representative parameter combo, grouped by starting option.
Returns a `Dict(option => Vector{Int})`.

Pathways that pause seeding — `:nothing` in *any* decision block, starting block included — are
excluded: the robustness analysis compares continuously-seeded pathways only.
"""
function group_pathways_by_starting_option(
    rs, option_names; sel_rcp, sel_dhw, sel_N_seed, sel_n_locations,
    seed_year_start, seed_years, pd_frequency, N_seed_weights
)
    max_time = size(rs.seed_log, :timesteps)
    combo_mask =
        rs.inputs.N_seed_CA .== sel_N_seed * N_seed_weights.N_seed_CA .&&
        rs.inputs.N_seed_TA .== sel_N_seed * N_seed_weights.N_seed_TA .&&
        rs.inputs.min_iv_locations .== sel_n_locations .&&
        rs.inputs.dhw_scenario .== sel_dhw .&&
        rs.inputs.RCP .== sel_rcp
    combo_idxs = findall(combo_mask)

    seqs = Dict(
        i => compact_pathway(
            rs.inputs.option_ts[i], seed_year_start, seed_years, pd_frequency, max_time
        )
        for i in combo_idxs
    )
    seeded_idxs = [i for i in combo_idxs if !any(==(:nothing), seqs[i])]

    return Dict(
        option => [i for i in seeded_idxs if seqs[i][1] == option]
        for option in option_names
    )
end

"""
    find_counterfactual(rs, sel_dhw, sel_rcp) -> Int

Column in `rs` of the no-seeding counterfactual (zero total seeds) for the given dhw + rcp.
pd_main.jl writes one such counterfactual per dhw_scenario.
"""
function find_counterfactual(rs, sel_dhw, sel_rcp)
    total_N_seed = rs.inputs.N_seed_TA .+ rs.inputs.N_seed_CA .+ rs.inputs.N_seed_CNA .+
                   rs.inputs.N_seed_SM .+ rs.inputs.N_seed_LM
    cf_idx = findfirst(
        (total_N_seed .== 0) .&
        (rs.inputs.dhw_scenario .== sel_dhw) .&
        (rs.inputs.RCP .== sel_rcp)
    )
    @assert !isnothing(cf_idx) "No counterfactual (zero total seeds) found in `rs` for dhw $(sel_dhw), RCP $(sel_rcp)."
    return cf_idx
end

# ── Aggregated per-pathway robustness (pathways analysis) ──────────────────────

"""
    pathway_robustness(opt_metric_mats, cf_metric_full, s; tail_fraction, seed_years)

Single scalar summarising pathway (scenario) `s`: the mean of the six tail ratios (bottom and
top of each of the three metrics, as returned by [`delta_tail_ratio`]). Metric 1
(`n_yrs_above`) uses `norm=seed_years` as its zero-guard, matching Figures A/B; metrics 2 and 3
use the default counterfactual-mean norm.
"""
function pathway_robustness(opt_metric_mats, cf_metric_full, s; tail_fraction, seed_years)
    b1, t1, _ = delta_tail_ratio(opt_metric_mats[1][:, s], cf_metric_full[1]; tail_fraction, norm=seed_years)
    b2, t2, _ = delta_tail_ratio(opt_metric_mats[2][:, s], cf_metric_full[2]; tail_fraction)
    b3, t3, _ = delta_tail_ratio(opt_metric_mats[3][:, s], cf_metric_full[3]; tail_fraction)
    return mean((b1, b2, b3, t1, t2, t3))
end

"""
    detect_param_sets(rs, sel_rcp; N_seed_weights) -> Vector{NamedTuple}

Unique `(N_seed, n_locations)` parameter sets among the *seeded* (non-counterfactual) scenarios
at `sel_rcp`, sorted ascending by `N_seed`. `N_seed` is the total seed budget reconstructed from
the per-taxa columns; counterfactual rows (zero budget) are dropped.
"""
function detect_param_sets(rs, sel_rcp; N_seed_weights)
    total_N_seed = rs.inputs.N_seed_TA .+ rs.inputs.N_seed_CA .+ rs.inputs.N_seed_CNA .+
                   rs.inputs.N_seed_SM .+ rs.inputs.N_seed_LM
    mask = (rs.inputs.RCP .== sel_rcp) .& (total_N_seed .> 0)
    combos = unique(
        [(N_seed=total_N_seed[i], n_locations=Int(rs.inputs.min_iv_locations[i]))
         for i in findall(mask)]
    )
    return sort(combos; by=c -> c.N_seed)
end

"""
    pathway_worst_dhw_robustness(rs, option_names, opt_metric_mats, ps; dhw_scenarios, sel_rcp,
        seed_year_start, seed_years, pd_frequency, N_seed_weights, tail_fraction)
        -> Dict(option => Dict(option_ts => Float64))

Worst-case robustness of every *individual* pathway of parameter set `ps`, grouped by starting
option. A pathway is identified by its `option_ts` code and its value is the minimum
[`pathway_robustness`] across the DHW scenarios in which it occurs, so the DHW reduction happens
per pathway and is independent of any starting-option aggregation.
"""
function pathway_worst_dhw_robustness(
    rs, option_names, opt_metric_mats, ps;
    dhw_scenarios, sel_rcp,
    seed_year_start, seed_years, pd_frequency, N_seed_weights, tail_fraction
)
    K = eltype(rs.inputs.option_ts)
    worst = Dict(option => Dict{K,Float64}() for option in option_names)

    for sel_dhw in dhw_scenarios
        pathways = group_pathways_by_starting_option(
            rs, option_names;
            sel_rcp, sel_dhw, sel_N_seed=ps.N_seed, sel_n_locations=ps.n_locations,
            seed_year_start, seed_years, pd_frequency, N_seed_weights
        )
        cf_idx = find_counterfactual(rs, sel_dhw, sel_rcp)
        cf_full = (
            opt_metric_mats[1][:, cf_idx],
            opt_metric_mats[2][:, cf_idx],
            opt_metric_mats[3][:, cf_idx]
        )

        for option in option_names
            d = worst[option]
            for s in pathways[option]
                key = rs.inputs.option_ts[s]
                rob = pathway_robustness(opt_metric_mats, cf_full, s; tail_fraction, seed_years)
                d[key] = haskey(d, key) ? min(d[key], rob) : rob
            end
        end
    end

    return worst
end

"""
    worst_dhw_robustness(rs, option_names, opt_metric_mats; dhw_scenarios, param_sets, sel_rcp,
        seed_year_start, seed_years, pd_frequency, N_seed_weights, tail_fraction) -> DataFrame

For each (starting option × parameter set), summarise the per-pathway worst-over-DHW robustness
([`pathway_worst_dhw_robustness`]) with its P10, median and P90 over pathways.

Returns a tidy table with columns `start_option, N_seed, n_locations, p10, median, p90`.
Options with no pathways in a parameter set are skipped.
"""
function worst_dhw_robustness(
    rs, option_names, opt_metric_mats;
    dhw_scenarios, param_sets, sel_rcp,
    seed_year_start, seed_years, pd_frequency, N_seed_weights, tail_fraction
)
    df = DataFrame(;
        start_option=String[], N_seed=Float64[], n_locations=Int[],
        p10=Float64[], median=Float64[], p90=Float64[]
    )

    for ps in param_sets
        worst = pathway_worst_dhw_robustness(
            rs, option_names, opt_metric_mats, ps;
            dhw_scenarios, sel_rcp,
            seed_year_start, seed_years, pd_frequency, N_seed_weights, tail_fraction
        )

        for option in option_names
            robs = collect(values(worst[option]))
            isempty(robs) && continue
            push!(df, (
                string(option), ps.N_seed, ps.n_locations,
                quantile(robs, 0.1), median(robs), quantile(robs, 0.9)
            ))
        end
    end

    return df
end

"""
    join_robustness_diversity(rob_df, pd_df) -> DataFrame

Pair each (starting option × parameter set) with a robustness value and a pathway-diversity
value. Robustness is `rob_df.median` (median of the per-pathway worst-over-DHW robustness from
[`worst_dhw_robustness`]), carried through with its P10/P90 spread over pathways; pathway
diversity is the *worst* (minimum over DHW scenarios) value from `pd_df` (an options-format
table: `option_name, pathway_diversity, N_seed, dhw_scenario, n_locations, rcp`).
Returns columns `start_option, N_seed, n_locations, robustness, robustness_p10, robustness_p90,
pathway_diversity, dominated?`.

`dominated?` flags Pareto-dominated options: within a parameter set (N_seed × n_locations), an
option is dominated when another option is at least as good in both dimensions (pathway diversity
and robustness) and strictly better in at least one.
"""
function join_robustness_diversity(rob_df, pd_df)
    worst = combine(
        groupby(pd_df, [:option_name, :N_seed, :n_locations])
    ) do sub
        (pathway_diversity=minimum(sub.pathway_diversity),)
    end
    worst.N_seed = Float64.(worst.N_seed)
    rename!(worst, :option_name => :start_option)

    joined = innerjoin(
        rob_df[:, [:start_option, :N_seed, :n_locations, :median, :p10, :p90]],
        worst; on=[:start_option, :N_seed, :n_locations]
    )
    rename!(joined, :median => :robustness, :p10 => :robustness_p10, :p90 => :robustness_p90)

    # Pareto dominance within each parameter set (N_seed × n_locations). An option is dominated if
    # some other option in the same panel is ≥ on both axes and strictly > on at least one.
    joined[!, "dominated?"] = falses(nrow(joined))
    for panel in groupby(joined, [:N_seed, :n_locations])
        x = panel.pathway_diversity
        y = panel.robustness
        for i in 1:nrow(panel)
            panel[i, "dominated?"] = any(
                j -> j != i && x[j] >= x[i] && y[j] >= y[i] && (x[j] > x[i] || y[j] > y[i]),
                1:nrow(panel)
            )
        end
    end
    return joined
end

# ── Rank agreement between robustness and pathway diversity ────────────────────

"""
    _kendall_counts(x, y) -> (nc, nd, tx, ty)

Count concordant (`nc`) and discordant (`nd`) pairs of `(x, y)`, plus pairs tied on `x` only
(`tx`) and on `y` only (`ty`). Pairs tied on both are ignored. These are the ingredients of the
Kendall τ-b statistic and, crucially, are *additive across strata* — summing them over panels
yields a blocked τ that only ever compares points within the same panel.
"""
function _kendall_counts(x::AbstractVector, y::AbstractVector)
    n = length(x)
    nc = nd = tx = ty = 0
    for i in 1:(n - 1), j in (i + 1):n
        dx = sign(x[i] - x[j])
        dy = sign(y[i] - y[j])
        if dx == 0 && dy == 0
            continue            # tied on both: contributes to neither τ-b denominator factor
        elseif dx == 0
            tx += 1             # tied on x only
        elseif dy == 0
            ty += 1             # tied on y only
        elseif dx == dy
            nc += 1
        else
            nd += 1
        end
    end
    return (nc, nd, tx, ty)
end

# τ-b from raw pair counts; NaN when one variable has no untied pairs (denominator = 0).
function _tau_b(nc, nd, tx, ty)
    denom = sqrt((nc + nd + ty) * (nc + nd + tx))
    return denom == 0 ? NaN : (nc - nd) / denom
end

"""
    kendall_tau_b(x, y) -> Float64

Kendall τ-b rank correlation between vectors `x` and `y`, with tie correction. Interpretable as
the (concordant − discordant) fraction of comparable pairs: +1 = identical ranking, −1 = exactly
reversed, 0 = no rank agreement. Returns `NaN` when either variable is all ties.
"""
kendall_tau_b(x::AbstractVector, y::AbstractVector) = _tau_b(_kendall_counts(x, y)...)

"""
    panel_tau_diversity_robustness(sub) -> Float64

Within-panel Kendall τ-b between pathway diversity and robustness, folding in the P10/median/P90
robustness spread as the mean of three separate τ-b's, one per robustness level.
"""
function panel_tau_diversity_robustness(sub)
    x = sub.pathway_diversity
    return mean((
        kendall_tau_b(x, sub.robustness_p10),
        kendall_tau_b(x, sub.robustness),
        kendall_tau_b(x, sub.robustness_p90),
    ))
end

"""
    kendall_tau_diversity_robustness(rob_div_df; group_cols=[:N_seed, :n_locations])
        -> (per_panel, stratified)

Rank agreement between pathway diversity and robustness across starting options, using the full
P10/median/P90 robustness spread (see [`panel_tau_diversity_robustness`]).

`per_panel` is a DataFrame (`group_cols`, `tau`, `n`) with a within-panel Kendall τ-b per
parameter set — the descriptive concordance shown on each facet of Figure D. `stratified` is a
single *blocked* τ-b: a blocked τ-b is computed per robustness level (concordant/discordant/tie
counts accumulated within each panel and summed across panels) and the three are averaged. Options
are only ever compared against options sharing the same `N_seed × n_locations` parameter set. This
holds the parameter-set gradient fixed and avoids the Simpson's-paradox risk of a naive pool.
"""
function kendall_tau_diversity_robustness(rob_div_df; group_cols=[:N_seed, :n_locations])
    per_panel = combine(groupby(rob_div_df, group_cols)) do sub
        (tau=panel_tau_diversity_robustness(sub), n=nrow(sub))
    end

    # Blocked τ-b: sum pair counts within panels, only comparing options in the same panel.
    _blocked_tau(yfun) = begin
        nc = nd = tx = ty = 0
        for sub in groupby(rob_div_df, group_cols)
            c = _kendall_counts(yfun(sub)...)
            nc += c[1]; nd += c[2]; tx += c[3]; ty += c[4]
        end
        _tau_b(nc, nd, tx, ty)
    end

    stratified = mean((
        _blocked_tau(s -> (s.pathway_diversity, s.robustness_p10)),
        _blocked_tau(s -> (s.pathway_diversity, s.robustness)),
        _blocked_tau(s -> (s.pathway_diversity, s.robustness_p90)),
    ))

    return (per_panel=per_panel, stratified=stratified)
end

# ── Top robustness pathways (pathways analysis) ────────────────────────────────

"""
    top_robustness_pathways(rs, option_names, opt_metric_mats; dhw_scenarios, param_sets,
        sel_rcp, seed_year_start, seed_years, pd_frequency, N_seed_weights, tail_fraction,
        n_top=10) -> DataFrame

Per (starting option × parameter set), rank the distinct per-block option sequences
([`compact_pathway`]) by their worst-over-DHW robustness — the minimum [`pathway_robustness`]
across the DHW scenarios in which the sequence occurs — and keep the top `n_top`.

Returns a table with columns `start_option, N_seed, n_locations, Top1 … TopN`, where each `TopK`
cell is the `" > "`-joined option sequence (empty string when fewer than `n_top` sequences exist).
"""
function top_robustness_pathways(
    rs, option_names, opt_metric_mats;
    dhw_scenarios, param_sets, sel_rcp,
    seed_year_start, seed_years, pd_frequency, N_seed_weights, tail_fraction, n_top::Int=10
)
    max_time = size(rs.seed_log, :timesteps)
    top_cols = [Symbol("Top$i") for i in 1:n_top]

    df = DataFrame(; start_option=String[], N_seed=Float64[], n_locations=Int[])
    for c in top_cols
        df[!, c] = String[]
    end

    for ps in param_sets
        worst = pathway_worst_dhw_robustness(
            rs, option_names, opt_metric_mats, ps;
            dhw_scenarios, sel_rcp,
            seed_year_start, seed_years, pd_frequency, N_seed_weights, tail_fraction
        )

        for option in option_names
            ranked = sort(collect(worst[option]); by=last, rev=true)
            tops = fill("", n_top)
            for (i, (option_ts, _)) in enumerate(ranked)
                i > n_top && break
                seq = compact_pathway(
                    option_ts, seed_year_start, seed_years, pd_frequency, max_time
                )
                tops[i] = join(string.(seq), " > ")
            end
            push!(df, (string(option), ps.N_seed, ps.n_locations, tops...))
        end
    end

    return df
end

# ── Probability-weighted tail statistics (pathways analysis) ───────────────────

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
        Mminus[p], Mplus[p], _ = delta_tail_ratio(
            intervention[:, p], counterfactual; tail_fraction=tail_frac, norm=reference
        )
    end

    return (
        top=pathway_tail_stats(Mplus, probs),
        bottom=pathway_tail_stats(Mminus, probs),
        n_eff=1 / sum(probs .^ 2)
    )
end

"""
    load_pathway_probabilities(prob_csv_path; sel_N_seed, sel_dhw, sel_n_locations, sel_rcp)

Read pathway adoption probabilities for the representative combo, returning a
`Dict(option_ts => probability)`.
"""
function load_pathway_probabilities(
    prob_csv_path; sel_N_seed, sel_dhw, sel_n_locations, sel_rcp
)
    prob_df = CSV.read(prob_csv_path, DataFrame)
    prob_df = prob_df[
        (prob_df.N_seed .== round(Int, sel_N_seed)) .&
        (prob_df.dhw_scenario .== sel_dhw) .&
        (prob_df.n_locations .== sel_n_locations) .&
        (prob_df.rcp .== sel_rcp), :]
    return Dict(row.option_ts => row.probability for row in eachrow(prob_df))
end

"""
    compute_weighted_tail_stats(option_names, option_pathways, weighted_metrics, prob_map, rs;
        tail_fraction)

Tidy per (starting option × metric × tail) table of probability-weighted tail statistics.
`weighted_metrics` is a vector of `(name, intervention_matrix, counterfactual_column, reference)`
tuples. `metric_idx` 1–6 matches `metric_labels` ordering (worst/best interleaved per metric).
"""
function compute_weighted_tail_stats(
    option_names, option_pathways, weighted_metrics, prob_map, rs; tail_fraction
)
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

    return weighted_tail_stats
end
