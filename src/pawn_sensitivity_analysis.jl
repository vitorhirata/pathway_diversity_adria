#=
PAWN sensitivity analysis on local tail metrics and GBR-scale metrics.

- factors are Sobol'-sampled from their distributions (`ADRIA.sample`) rather than crossed
  over hand-picked levels; only the bounds/levels of each distribution are specified here
- the seeding option is crossed manually over the Sobol set (it is not a model-spec factor)
- N_seed split by fixed group weights (N_seed_weights), as in
  `static_options.jl`, instead of an equal 1/5 split
- seeding_devices_per_m2 and a_adapt are varied factors
- a counterfactual (no-seeding) scenario is added per dhw member
- the PAWN output metrics are, at two scales:
  * local tail metrics, counterfactual-relative, as in
    `static_options.jl`: for each of cumulative cover, years above 20%
    cover, and cumulative evenness, the per-reef delta (option − counterfactual) is reduced
    to the mean of the bottom/top `tail_number` reefs, normalized (worst/best reefs)
  * GBR-scale metrics: total absolute cover, relative shelter volume, relative juveniles and
    coral evenness, over all locations, reduced to their time mean. Reported twice: raw, and
    as the fractional change against the counterfactual. Both drop the no-seeding scenarios
    before PAWN, so the raw version measures how the factors move the absolute state and the
    relative version how they move the benefit of seeding

Varied factors:
- dhw_scenario
- min_iv_locations
- N_seed (total, split by fixed weights)
- Seed weights (option)
- seeding_devices_per_m2
- a_adapt
- RCP
=#

include("src/common.jl")
using GeoMakie, GraphMakie, CairoMakie
using ADRIA.Distributions: LogUniform, TriangularDist  # quantile comes from Statistics (extended)

# ── Constants ─────────────────────────────────────────────────────────────────

seed_years = 20
tail_number = 150
rcps = ["26", "45", "70"]
n_sobol = 512  # must be a power of 2 (Sobol' sampler requirement)

# Distribution bounds/levels of the sampled factors. dhw_scenario is left at its default
# category tuple, i.e. every DHW member in the domain.
N_seed_bounds = (1e6, 1e8)      # total corals per deployment, sampled log-triangular (mode below)
N_seed_peak = 3e6               # mode of the log-triangular budget: dense around realistic sizes
min_loc_bounds = (50.0, 1000.0) # number of intervention locations, sampled log-uniform
device_levels = (2.0, 5.0, 10.0)
a_adapt_bounds = (2.0, 10.0, 0.5)  # (lower, upper, step)
N_seed_weights = (
    N_seed_TA=0.15, N_seed_CA=0.5, N_seed_CNA=0.0, N_seed_SM=0.35, N_seed_LM=0.0
)

dom = ADRIA.load_domain(
    pd_config["domain_path"], rcps[1];
    calib_params_fn=pd_config["coral_param_path"],
    # timeframe: seed_years + 2 (start seeding), 5 (extra years)
    timeframe=(2022, 2022 + seed_years + 2 + 5)
)
fix_common_parameters!(dom)

ms = ADRIA.model_spec(dom)

# seeding_devices_per_m2, a_adapt, min_iv_locations, N_seed_CA and dhw_scenario are
# intentionally NOT fixed here: they are the sampled factors.
ADRIA.fix_factor!(dom; seed_years=seed_years)

# Only N_seed_CA is sampled; it carries the *total* budget and is split into groups below.
# Its lower bound must stay > 0, otherwise `adjust_samples` treats the row as unseeded and
# zeroes a_adapt and every seed_* column.
ADRIA.fix_factor!(dom; N_seed_TA=0, N_seed_CNA=0, N_seed_SM=0, N_seed_LM=0)

# Seed criteria weights are set per option below, so they must not be design dimensions
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "SeedCriteriaWeights").fieldname)

ADRIA.set_factor_bounds!(
    dom;
    N_seed_CA=(N_seed_bounds..., 50_000.0),
    min_iv_locations=min_loc_bounds,
    seeding_devices_per_m2=device_levels,
    a_adapt=a_adapt_bounds
)

# Beyond the five sampled factors, the remaining non-constant factors are the inert
# fog_*/mc_*/shade_*/reactive_*/mcb_deployment_freq group, which `adjust_samples` zeroes out
# given fogging=0, N_mc_settlers=0, SRM=0 and a periodic seed_strategy.
ms = ADRIA.model_spec(dom)
@info "Sobol design dimensions: $(ms[.!ms.is_constant, :fieldname])"

options = ADRIA.analysis.option_seed_preference(; include_weights=true)
n_opts = nrow(options)

# ── Build scenario table ──────────────────────────────────────────────────────

# Every option is paired with the identical Sobol set
scens = repeat(ADRIA.sample(dom, n_sobol); inner=n_opts)
scens[!, :option] = repeat(1:n_opts, n_sobol)

# Log-uniform number of intervention locations. `min_iv_locations` is sampled uniformly over a
# regular integer grid, so `u_loc` is uniform on [0, 1] and `quantile(LogUniform(...))` is a
# valid inverse-CDF transform. Concentrates samples toward smaller footprints, where adding a
# location changes per-location seed density (and its marginal effect) fastest.
u_loc = (scens.min_iv_locations .- min_loc_bounds[1]) ./ (min_loc_bounds[2] - min_loc_bounds[1])
scens[!, :min_iv_locations] = round.(Int64, quantile.(LogUniform(min_loc_bounds...), u_loc))

# Log-triangular total seeding budget. The sampled value is uniform over a regular grid, so `u`
# is uniform on [0, 1] and applying a triangular inverse-CDF in log10-space is a valid transform.
# The mode sits at `N_seed_peak`, keeping samples dense around realistic deployment sizes while
# the linear upper tail still reaches the 1e8 feasibility limit (~⅓ of samples land above 1e7 —
# the price of retaining 1e8 samples for the rest of the analysis).
u = (scens.N_seed_CA .- N_seed_bounds[1]) ./ (N_seed_bounds[2] - N_seed_bounds[1])
log_budget = TriangularDist(log10(N_seed_bounds[1]), log10(N_seed_bounds[2]), log10(N_seed_peak))
N_seed_total = 10.0 .^ quantile.(log_budget, u)

for (group, weight) in pairs(N_seed_weights)
    scens[!, group] = N_seed_total .* weight
end
# NOTE: this bookkeeping column must not contain "N_seed" in its name — `run_model` selects
# the seeding volume by substring (`ADRIA.jl/src/scenario.jl:909`), so an extra "N_seed*"
# column silently joins the per-group seed volumes and breaks the run.
scens[!, :total_deployed_corals] = N_seed_total

# Assigned by name: both `options` and the model spec derive their criteria from
# `fieldnames(ADRIA.SeedCriteriaWeights)`, so this stays correct if criteria are added or
# removed (positional indexing into `options` silently breaks when they are).
criteria = ADRIA.component_params(ms, "SeedCriteriaWeights").fieldname
for (row, scen) in enumerate(eachrow(scens))
    scens[row, criteria] = collect(options[scen.option, criteria])
end

# Counterfactual (no seeding) — one per dhw member. dhw_scenario is the only non-intervention
# factor that varies, so this is a complete set (RCP is crossed by `run_scenarios`). The other
# seeding parameters (seeding_devices_per_m2, a_adapt, option, min_iv_locations) are left at
# the first sampled value since they are irrelevant with no seeding; these rows are excluded
# from PAWN.
dhw_members = sort(unique(scens.dhw_scenario))
cf = repeat(scens[1:1, :], length(dhw_members))
cf[!, :dhw_scenario] = dhw_members
cf[!, :option] .= 0
cf[!, :total_deployed_corals] .= 0.0
for group in keys(N_seed_weights)
    cf[!, group] .= 0.0
end

scens = vcat(scens, cf)

# Run one ResultSet per RCP, adding RCP as a factor column, then combine
rs = ADRIA.run_scenarios(dom, scens, rcps)

# Load scenario
# path = "Outputs/"
# rs = ADRIA.load_results(path)

# ── Identify counterfactual vs intervention scenarios ─────────────────────────

total_N_seed = rs.inputs.N_seed_TA .+ rs.inputs.N_seed_CA .+ rs.inputs.N_seed_CNA .+
               rs.inputs.N_seed_SM .+ rs.inputs.N_seed_LM
cf_idxs = findall(total_N_seed .== 0)
intervention_idxs = findall(total_N_seed .> 0)

# ── Local metrics over full timeframe (computed once for all scenarios) ────────

n_locs = size(rs.seed_log, :locations)
n_scens_rs = nrow(rs.inputs)

m_tac = Array(ADRIA.metrics.total_absolute_cover(rs)) .* 1e-6  # (timesteps, locations, scenarios)
loc_hab_area_km2 = rs.loc_area .* rs.loc_max_coral_cover .* 1e-6
fd_data = ADRIA.metrics.coral_evenness(rs)                     # (timesteps, locations, scenarios)

n_yrs_above = ADRIA.ZeroDataCube((:locations, :scenarios), (n_locs, n_scens_rs); T=Int32)
for s in 1:n_scens_rs, l in 1:n_locs
    thr = 0.20 * loc_hab_area_km2[l]
    n_yrs_above.data[l, s] = Int32(count(m_tac[:, l, s] .>= thr))
end

cum_tac = dropdims(sum(m_tac; dims=1); dims=1)   # (locations, scenarios)
cum_fd = dropdims(sum(fd_data; dims=:timesteps); dims=:timesteps)  # (locations, scenarios)

# ── Counterfactual-relative tail metric ───────────────────────────────────────

"""
Select the bottom and top `tail_number` reefs by delta (opt − cf), then compute
mean(delta[group]) / norm for each group. If `norm == 0`, defaults to mean(cf).
Returns (bottom_ratio, top_ratio).
"""
function delta_tail_ratio(
    opt::AbstractVector, cf::AbstractVector; tail_number::Int=150, norm::Float64=0.0
)
    delta = opt .- cf
    order = sortperm(delta)
    bot = order[1:tail_number]
    top = order[(end - tail_number + 1):end]

    if iszero(norm)
        norm = mean(cf)
    end

    bot_ratio = mean(delta[bot]) / norm
    top_ratio = mean(delta[top]) / norm
    return bot_ratio, top_ratio
end

# Match each intervention scenario to its counterfactual by (dhw_scenario, RCP)
cf_lookup = Dict(
    (rs.inputs.dhw_scenario[c], rs.inputs.RCP[c]) => c for c in cf_idxs
)

# Per-scenario tail scalars (interventions filled, counterfactual entries left as NaN)
y_nyrs_bot = fill(NaN, n_scens_rs); y_nyrs_top = fill(NaN, n_scens_rs)
y_tac_bot = fill(NaN, n_scens_rs);  y_tac_top = fill(NaN, n_scens_rs)
y_fd_bot = fill(NaN, n_scens_rs);   y_fd_top = fill(NaN, n_scens_rs)

for s in intervention_idxs
    cf_s = cf_lookup[(rs.inputs.dhw_scenario[s], rs.inputs.RCP[s])]

    opt_nyrs = Float64.(n_yrs_above.data[:, s])
    cf_nyrs = Float64.(n_yrs_above.data[:, cf_s])
    opt_tac = cum_tac[:, s]
    cf_tac = cum_tac[:, cf_s]
    opt_fd = cum_fd.data[:, s]
    cf_fd = cum_fd.data[:, cf_s]

    y_nyrs_bot[s], y_nyrs_top[s] =
        delta_tail_ratio(opt_nyrs, cf_nyrs; tail_number=tail_number, norm=1.0)
    y_tac_bot[s], y_tac_top[s] =
        delta_tail_ratio(opt_tac, cf_tac; tail_number=tail_number, norm=1.0)
    y_fd_bot[s], y_fd_top[s] =
        delta_tail_ratio(opt_fd, cf_fd; tail_number=tail_number, norm=1.0)
end

# ── GBR-scale metrics, raw and counterfactual-relative ────────────────────────

"""
Fractional change of a per-scenario scalar against its matched counterfactual,
`(option − counterfactual) / counterfactual`. Counterfactual entries are left as NaN and are
dropped before PAWN.
"""
function cf_relative_change(scen_metric::AbstractVector)::Vector{Float64}
    out = fill(NaN, length(scen_metric))
    for s in intervention_idxs
        cf_s = cf_lookup[(rs.inputs.dhw_scenario[s], rs.inputs.RCP[s])]
        out[s] = (scen_metric[s] - scen_metric[cf_s])
    end
    return out
end

# Scenario-level metrics over all locations, reduced to their mean over time. These are the
# raw PAWN outputs; the no-seeding scenarios are dropped downstream by `intervention_idxs`.
y_s_tac_raw = vec(mean(ADRIA.metrics.scenario_total_cover(rs); dims=:timesteps))
y_s_rsv_raw = vec(mean(ADRIA.metrics.scenario_rsv(rs); dims=:timesteps))
y_s_juves_raw = vec(mean(ADRIA.metrics.scenario_relative_juveniles(rs); dims=:timesteps))
y_s_even_raw = vec(mean(ADRIA.metrics.scenario_evenness(rs); dims=:timesteps))

# GBR-scale years above 20% cover: sum the per-reef counts over locations (kept in YAXArray
# form, then reduced to a per-scenario vector).
n_yrs_above_gbr = dropdims(sum(n_yrs_above; dims=:locations); dims=:locations)
y_s_nyrs_raw = Float64.(n_yrs_above_gbr.data)

y_s_tac = cf_relative_change(y_s_tac_raw)
y_s_rsv = cf_relative_change(y_s_rsv_raw)
y_s_juves = cf_relative_change(y_s_juves_raw)
y_s_even = cf_relative_change(y_s_even_raw)
y_s_nyrs = cf_relative_change(y_s_nyrs_raw)

# ── PAWN sensitivity per metric × tail ────────────────────────────────────────

fig_opts = Dict(:size => (1500, 1000))
opts = Dict(
    :factors => [
        :dhw_scenario, :dhw_mean, :RCP, :option, :min_iv_locations, :total_deployed_corals,
        :seeding_devices_per_m2, :a_adapt
    ],
    :by => :none,
    :ytick_labels => [
        "Climate Model",
        "Mean DHW",
        "RCP",
        "Option",
        "Number of locations",
        "Number of corals deployed",
        "Seeding devices per m²",
        "Assisted adaptation"
    ]
)
axis_opts = Dict(
    :title => "",
    :titlesize => 22,
    :xticklabelsize => 20,
    :yticklabelsize => 20
)

# Each figure groups one metric type into a 2×2 panel of PAWN indices sharing a single
# colorbar. Panels are laid out row-major: [(1,1), (1,2), (2,1), (2,2)].
figure_specs = [
    (
        suptitle="Cumulative cover",
        filename="pawn_cumulative_cover",
        panels=[
            (y_s_tac_raw, "GBR (raw)"),
            (y_s_tac, "GBR (vs counterfactual)"),
            (y_tac_top, "Best reefs (vs counterfactual)"),
            (y_tac_bot, "Worst reefs (vs counterfactual)")
        ]
    ),
    (
        suptitle="Cumulative evenness",
        filename="pawn_cumulative_evenness",
        panels=[
            (y_s_even_raw, "GBR (raw)"),
            (y_s_even, "GBR (vs counterfactual)"),
            (y_fd_top, "Best reefs (vs counterfactual)"),
            (y_fd_bot, "Worst reefs (vs counterfactual)")
        ]
    ),
    (
        suptitle="Years above 20% cover",
        filename="pawn_years_above",
        panels=[
            (y_s_nyrs_raw, "GBR (raw)"),
            (y_s_nyrs, "GBR (vs counterfactual)"),
            (y_nyrs_top, "Best reefs (vs counterfactual)"),
            (y_nyrs_bot, "Worst reefs (vs counterfactual)")
        ]
    ),
    (
        suptitle="Shelter volume and relative juveniles",
        filename="pawn_shelter_juveniles",
        panels=[
            (y_s_rsv_raw, "Relative shelter volume (GBR, raw)"),
            (y_s_rsv, "Relative shelter volume (GBR, vs counterfactual)"),
            (y_s_juves_raw, "Relative juveniles (GBR, raw)"),
            (y_s_juves, "Relative juveniles (GBR, vs counterfactual)")
        ]
    )
]

"""
Mean DHW of each ensemble member, per RCP.

`dhw_scenario` is an arbitrary member ID, so PAWN's quantile slices over it group unrelated
members together and dilute the index. Slicing on mean DHW instead orders the members by
the heat stress they actually impose. Both are reported.

Computed from the domain rather than `rs.dhw_stats`: `setup_result_store!` builds the per-RCP
DHW summaries before `run_scenarios` switches RCPs, so every entry there holds the stats of
whichever RCP the domain was loaded with.
"""
dhw_mean_by_rcp = Dict{Int64,Vector{Float64}}()
for rcp in rcps
    ADRIA.switch_RCPs!(dom, rcp)
    # dims: (timesteps, locations, members)
    dhw_mean_by_rcp[parse(Int64, rcp)] = vec(mean(dom.dhw_scens.data; dims=(1, 2)))
end

X = rs.inputs[intervention_idxs, :]
X[!, :dhw_mean] = [
    dhw_mean_by_rcp[Int64(scen.RCP)][Int64(scen.dhw_scenario)] for scen in eachrow(X)
]

for spec in figure_specs
    f = Figure(; fig_opts...)
    for (i, (y, ptitle)) in enumerate(spec.panels)
        row = (i - 1) ÷ 2 + 1
        col = (i - 1) % 2 + 1
        si = ADRIA.sensitivity.pawn(X, y[intervention_idxs])
        g = f[row, col] = GridLayout()
        ax_opts = copy(axis_opts)
        ax_opts[:title] = ptitle
        # Factor labels are identical across panels; keep them only on the left column.
        if col == 2
            ax_opts[:yticklabelsvisible] = false
        end
        ADRIA.viz.pawn!(g, si; opts=opts, axis_opts=ax_opts)
        Label(
            g[1, 1, TopLeft()], "$('a' + i - 1))";
            fontsize=26,
            font=:bold,
            halign=:right,
            padding=(0, 8, 8, 0)
        )
    end
    Colorbar(
        f[1:2, 3];
        colormap=:viridis,
        colorrange=(0, 1),
        label="Relative Sensitivity",
        labelsize=22,
        ticklabelsize=20,
        height=Relative(0.65)
    )
    Label(f[0, 1:2], spec.suptitle; fontsize=28, font=:bold)
    save(joinpath(pd_config["plot_output_path"], "$(spec.filename).png"), f)
end
