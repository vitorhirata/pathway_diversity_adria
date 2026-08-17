#=
Data processing for the static-option detail plots (maps, seeding-frequency histograms,
per-option time-series).

These plots were originally written for a single-parameter run laid out as fixed blocks of
`n_scens_per_dhw` rows. In the merged `static_options.jl` they run on top of the robustness
parameter-sweep scenario table (RCP × N_seed × min_iv_locations × DHW × 7 rows), so scenarios are
selected by matching `rs.inputs` columns and the `:option` identifier (-1 = counterfactual,
0 = unguided, 1:5 = guided option), the same idiom as `cvar_aggregation` in
`data_processing/robustness.jl`.

This file is *included by* the main script and assumes `include("src/common.jl")` has already run
(it uses `ADRIA`, `Statistics`, `DataFrames`). No plotting lives here — see
`visualization/static_options.jl`.
=#

"""
    select_static_scenarios(rs, N_seed_total, min_iv, rcp;
                            N_seed_weights, dhw_scenarios, scenario_names)

Resolve the `rs` scenario indices for one static-option parameter set. Returns a NamedTuple `sel`:
- `opt_cols`  : `n_dhw × n_options` matrix of `rs` scenario indices, option order = `scenario_names`
  (last option is `:unguided`, `option == 0`).
- `cf_cols`   : per-DHW counterfactual (`option == -1`) scenario index.
- `intervention_scens` : flat vector of the option columns, ordered so `iv_col(d, o)` indexes it.
- `iv_col(d, o)`  : column of option `o` under DHW block `d` in intervention-only arrays.
- `block_cols(d)` : the 7 columns of DHW block `d` in the order `[options..., counterfactual]`,
  matching the per-DHW block layout the time-series plot expects.
"""
function select_static_scenarios(
    rs, N_seed_total, min_iv, rcp;
    N_seed_weights, dhw_scenarios, scenario_names
)
    n_options = length(scenario_names)
    n_dhw = length(dhw_scenarios)
    rcp_val = parse(Float64, rcp)

    # Common filter for this parameter set; the counterfactual carries no N_seed so it is matched
    # by `param_mask` alone, while seeded rows additionally match the seed budget.
    param_mask = (rs.inputs.RCP .== rcp_val) .& (rs.inputs.min_iv_locations .== min_iv)
    seeded_mask = param_mask .& (rs.inputs.N_seed_CA .== N_seed_total * N_seed_weights.N_seed_CA)

    opt_cols = Matrix{Int}(undef, n_dhw, n_options)
    cf_cols = Vector{Int}(undef, n_dhw)
    for (d, dhw) in enumerate(dhw_scenarios)
        dhw_mask = rs.inputs.dhw_scenario .== dhw
        cf = findfirst(param_mask .& dhw_mask .& (rs.inputs.option .== -1))
        @assert !isnothing(cf) "No counterfactual for rcp $(rcp), min_iv $(min_iv), dhw $(dhw)."
        cf_cols[d] = cf
        for o in 1:n_options
            # option id: 1:5 for guided options, 0 for the unguided (last) column.
            opt_id = o < n_options ? o : 0
            s = findfirst(seeded_mask .& dhw_mask .& (rs.inputs.option .== opt_id))
            @assert !isnothing(s) "No option $(scenario_names[o]) for rcp $(rcp), N_seed $(N_seed_total), min_iv $(min_iv), dhw $(dhw)."
            opt_cols[d, o] = s
        end
    end

    intervention_scens = vcat([opt_cols[d, :] for d in 1:n_dhw]...)
    iv_col = (d, o) -> (d - 1) * n_options + o
    block_cols = d -> vcat(opt_cols[d, :], cf_cols[d])

    return (;
        opt_cols, cf_cols, intervention_scens, iv_col, block_cols,
        n_dhw, n_options, dhw_scenarios, scenario_names
    )
end

"""
    seeding_derived(rs, sel) -> (; seed_per_reef_per_ts_scen, total_seeds,
                                   always_seeded, never_seeded, active_mask, selected_locations)

Seeding quantities over the intervention scenarios of parameter set `sel`, scoped to the seeding
window. `seed_per_reef_per_ts_scen` and `total_seeds` are ordered so column `sel.iv_col(d, o)`
holds option `o` under DHW block `d`. `always/never_seeded` hold across every timestep and every
intervention scenario of the set.
"""
function seeding_derived(rs, sel)
    seed_start = Int(rs.inputs.seed_year_start[1])
    n_seed_years = Int(rs.inputs.seed_years[1])
    seed_ts = seed_start:(seed_start + n_seed_years - 1)

    seed_per_reef_per_ts_scen = dropdims(
        sum(rs.seed_log[timesteps=seed_ts, scenarios=sel.intervention_scens]; dims=:coral_id);
        dims=:coral_id
    )
    total_seeds = Array(
        dropdims(
            sum(rs.seed_log[scenarios=sel.intervention_scens]; dims=(:timesteps, :coral_id));
            dims=(:timesteps, :coral_id)
        )
    )  # (n_locs, n_dhw * n_options), column sel.iv_col(d, o)

    # always/never seeded: must hold across every timestep AND every intervention scenario
    always_seeded = vec(all(seed_per_reef_per_ts_scen.data .> 0; dims=(1, 3)))
    never_seeded = vec(all(seed_per_reef_per_ts_scen.data .== 0; dims=(1, 3)))
    active_mask = .!never_seeded
    selected_locations = findall(.!(always_seeded .| never_seeded))

    return (;
        seed_per_reef_per_ts_scen, total_seeds,
        always_seeded, never_seeded, active_mask, selected_locations
    )
end

"""
    static_performance_metrics(rs, sel) -> (; n_yrs_diff, cum_tac_diff)

Per-reef performance deltas vs the same-DHW counterfactual for parameter set `sel`. Both matrices
are `(n_locs, n_dhw * n_options)`, ordered by `sel.iv_col(d, o)`:
- `n_yrs_diff`   : Δ count of timesteps with cover ≥ 20% habitable area (option − counterfactual).
- `cum_tac_diff` : Δ cumulative absolute cover (km²·years, option − counterfactual).
"""
function static_performance_metrics(rs, sel)
    m_tac = Array(ADRIA.metrics.total_absolute_cover(rs)) .* 1e-6  # km²
    loc_hab_area_km2 = rs.loc_area .* rs.loc_max_coral_cover .* 1e-6
    n_locs = size(m_tac, 2)

    # Per-scenario per-reef metrics
    nyrs(s) = Float64[count(m_tac[:, l, s] .>= 0.20 * loc_hab_area_km2[l]) for l in 1:n_locs]
    cumtac(s) = vec(sum(m_tac[:, :, s]; dims=1))  # (n_locs,), km²·years

    n_cols = sel.n_dhw * sel.n_options
    n_yrs_diff = Matrix{Float64}(undef, n_locs, n_cols)
    cum_tac_diff = Matrix{Float64}(undef, n_locs, n_cols)
    for d in 1:sel.n_dhw
        cf_nyrs = nyrs(sel.cf_cols[d])
        cf_cumtac = cumtac(sel.cf_cols[d])
        for o in 1:sel.n_options
            col = sel.iv_col(d, o)
            n_yrs_diff[:, col] .= nyrs(sel.opt_cols[d, o]) .- cf_nyrs
            cum_tac_diff[:, col] .= cumtac(sel.opt_cols[d, o]) .- cf_cumtac
        end
    end

    return (; n_yrs_diff, cum_tac_diff)
end

"""
    scenario_timeseries_metrics(rs, selected_locations) -> Dict{String, YAXArray}

Scenario-level time-series metrics restricted to `selected_locations`, keyed by display name.
Columns are aligned with `rs`'s scenario dimension (index with `sel.block_cols(d)` to get a
DHW block). `scenario_relative_juveniles` ignores the `locations` kwarg, so it is pre-sliced.
"""
function scenario_timeseries_metrics(rs, selected_locations)
    s_tac = ADRIA.metrics.scenario_total_cover(rs; locations=selected_locations)
    s_rsv = ADRIA.metrics.scenario_rsv(rs; locations=selected_locations)
    s_even = ADRIA.metrics.scenario_evenness(rs; locations=selected_locations)

    _aj = ADRIA.metrics.absolute_juveniles(rs)
    _k_area = ADRIA.loc_k_area(rs)[selected_locations]
    s_juves = ADRIA.metrics.scenario_relative_juveniles(
        _aj[locations=selected_locations].data, _k_area
    )

    return Dict(
        "Total absolute cover" => s_tac,
        "Relative Shelter Volume" => s_rsv,
        "Relative Juveniles" => s_juves,
        "Coral Evenness" => s_even
    )
end
