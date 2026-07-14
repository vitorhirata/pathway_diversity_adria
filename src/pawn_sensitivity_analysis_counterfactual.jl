#=
PAWN sensitivity analysis on local, counterfactual-relative tail metrics.

Same PAWN scaffolding as `pawn_sensitivity_analysis.jl`, but:
- dhw_scenarios restricted to [2, 5, 7, 9, 11]
- N_seed split by fixed group weights (N_seed_weights), as in
  `robustness_analysis_static_options.jl`, instead of an equal 1/5 split
- seeding_devices_per_m2 and a_adapt are varied factors ([1, 5, 10])
- a counterfactual (no-seeding) scenario is added per dhw
- the PAWN output metric is the local, counterfactual-relative tail metric used in
  `robustness_analysis_static_options.jl`: for each of cumulative cover, years above 20%
  cover, and cumulative evenness, the per-reef delta (option − counterfactual) is reduced
  to the mean of the bottom/top `tail_fraction` of reefs, normalized (worst/best reefs).

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

# ── Constants ─────────────────────────────────────────────────────────────────

seed_years = 20
tail_fraction = 0.05
rcps = ["26", "45", "70"]
dhw_scenarios = [5, 7, 9, 11]

N_seeds = [1e6, 1e7, 1e8]
min_locations = [50, 200, 500, 1000]
seeding_devices = [2, 5, 10]
a_adapts = [2, 5, 10]
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

# seeding_devices_per_m2 and a_adapt are intentionally NOT fixed here: they are varied.
ADRIA.fix_factor!(dom; seed_years=seed_years)

options = ADRIA.analysis.option_seed_preference(; include_weights=true)

# ── Build scenario table ──────────────────────────────────────────────────────

n_intervention = length(N_seeds) * nrow(options) * length(dhw_scenarios) *
                 length(min_locations) * length(seeding_devices) * length(a_adapts)
n_cf = length(dhw_scenarios)  # one counterfactual per dhw
n_scens = n_intervention + n_cf

scens_base = repeat(ADRIA.sample(dom, 2)[1:1, :], n_scens)
scens_base[!, :option] = zeros(Int, n_scens)

row = 1
for N_seed in N_seeds, (opt_idx, option) in enumerate(eachrow(options)),
    dhw_scenario in dhw_scenarios, min_location in min_locations,
    device in seeding_devices, a_adapt in a_adapts

    scens_base[row, :N_seed_TA] = N_seed * N_seed_weights.N_seed_TA
    scens_base[row, :N_seed_CA] = N_seed * N_seed_weights.N_seed_CA
    scens_base[row, :N_seed_CNA] = N_seed * N_seed_weights.N_seed_CNA
    scens_base[row, :N_seed_SM] = N_seed * N_seed_weights.N_seed_SM
    scens_base[row, :N_seed_LM] = N_seed * N_seed_weights.N_seed_LM
    scens_base[row, :seed_heat_stress] = option[2]
    scens_base[row, :seed_in_connectivity] = option[3]
    scens_base[row, :seed_out_connectivity] = option[4]
    scens_base[row, :seed_depth] = option[5]
    scens_base[row, :seed_coral_cover] = option[6]
    scens_base[row, :seed_cluster_diversity] = option[7]
    scens_base[row, :seed_geographic_separation] = option[8]
    scens_base[row, :seed_coral_diversity] = option[9]
    scens_base[row, :dhw_scenario] = dhw_scenario
    scens_base[row, :min_iv_locations] = min_location
    scens_base[row, :seeding_devices_per_m2] = device
    scens_base[row, :a_adapt] = a_adapt
    scens_base[row, :option] = opt_idx
    row += 1
end

# Counterfactual (no seeding) — one per dhw. Other seeding parameters
# (seeding_devices_per_m2, a_adapt, option, min_iv_locations) are left at the base
# sample value since they are irrelevant with no seeding; these rows are excluded from PAWN.
for dhw_scenario in dhw_scenarios
    scens_base[row, :dhw_scenario] = dhw_scenario
    scens_base[row, :N_seed_TA] = 0
    scens_base[row, :N_seed_CA] = 0
    scens_base[row, :N_seed_CNA] = 0
    scens_base[row, :N_seed_SM] = 0
    scens_base[row, :N_seed_LM] = 0
    scens_base[row, :option] = 0
    row += 1
end

# Run one ResultSet per RCP, adding RCP as a factor column, then combine
rs = ADRIA.run_scenarios(dom, scens_base, rcps)

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

# Reefs seeded in at least one intervention scenario (across all options, dhw, RCPs)
seed_start = Int(rs.inputs.seed_year_start[1])
n_seed_years = Int(rs.inputs.seed_years[1])
seed_ts = seed_start:(seed_start + n_seed_years - 1)
seed_per_reef_per_ts_scen = dropdims(
    sum(rs.seed_log[timesteps=seed_ts, scenarios=intervention_idxs]; dims=:coral_id);
    dims=:coral_id
)
never_seeded = vec(all(seed_per_reef_per_ts_scen.data .== 0; dims=(1, 3)))
active_mask = .!never_seeded

# ── Counterfactual-relative tail metric ───────────────────────────────────────

"""
Select bottom and top `tail_fraction` of reefs by delta (opt − cf), then compute
mean(delta[group]) / norm for each group. If `norm == 0`, defaults to mean(cf).
Returns (bottom_ratio, top_ratio).
"""
function delta_tail_ratio(
    opt::AbstractVector, cf::AbstractVector; tail_fraction::Float64=0.05, norm::Float64=0.0
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

    opt_nyrs = Float64.(n_yrs_above.data[active_mask, s])
    cf_nyrs = Float64.(n_yrs_above.data[active_mask, cf_s])
    opt_tac = cum_tac[active_mask, s]
    cf_tac = cum_tac[active_mask, cf_s]
    opt_fd = cum_fd.data[active_mask, s]
    cf_fd = cum_fd.data[active_mask, cf_s]

    y_nyrs_bot[s], y_nyrs_top[s] =
        delta_tail_ratio(opt_nyrs, cf_nyrs; tail_fraction=tail_fraction, norm=float(seed_years))
    y_tac_bot[s], y_tac_top[s] =
        delta_tail_ratio(opt_tac, cf_tac; tail_fraction=tail_fraction)
    y_fd_bot[s], y_fd_top[s] =
        delta_tail_ratio(opt_fd, cf_fd; tail_fraction=tail_fraction)
end

# ── PAWN sensitivity per metric × tail ────────────────────────────────────────

fig_opts = Dict(:size => (800, 400))
opts = Dict(
    :factors => [
        :dhw_scenario, :RCP, :option, :min_iv_locations, :N_seed_TA,
        :seeding_devices_per_m2, :a_adapt
    ],
    :by => :none,
    :ytick_labels => [
        "Climate Model",
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

metric_ys = [
    y_tac_bot, y_tac_top,
    y_nyrs_bot, y_nyrs_top,
    y_fd_bot, y_fd_top
]
titles = [
    "Cumulative cover (worst reefs)", "Cumulative cover (best reefs)",
    "Years above 20% cover (worst reefs)", "Years above 20% cover (best reefs)",
    "Cumulative evenness (worst reefs)", "Cumulative evenness (best reefs)"
]
filenames = [
    "pawn_cf_cumulative_cover_worst", "pawn_cf_cumulative_cover_best",
    "pawn_cf_years_above_worst", "pawn_cf_years_above_best",
    "pawn_cf_cumulative_evenness_worst", "pawn_cf_cumulative_evenness_best"
]

X = rs.inputs[intervention_idxs, :]
for (idx, y) in enumerate(metric_ys)
    si = ADRIA.sensitivity.pawn(X, y[intervention_idxs])
    axis_opts[:title] = titles[idx]
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.pawn!(g, si; opts=opts, axis_opts=axis_opts)
    Colorbar(
        g[1, 2];
        colormap=:viridis,
        colorrange=(0, 1),
        label="Relative Sensitivity",
        labelsize=22,
        ticklabelsize=20,
        height=Relative(0.65)
    )
    save(joinpath(pd_config["plot_output_path"], "$(filenames[idx]).png"), f)
end
