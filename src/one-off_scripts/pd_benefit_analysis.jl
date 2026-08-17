include("src/common.jl")
using WGLMakie, GraphMakie, CairoMakie

path = ""
rs = ADRIA.load_results(path)

# ── User-configurable variables ───────────────────────────────────────────────

previous_pathway = [:heat_stress, :heat_stress, :heat_stress]

dhw_scenario = 11
N_seed = 1e7
min_iv_locations = 200
RCP = "45"

N_seed_weights = (
    N_seed_TA=0.15, N_seed_CA=0.5, N_seed_CNA=0.0, N_seed_SM=0.35, N_seed_LM=0.0
)

# ── Decode all option_ts ──────────────────────────────────────────────────────

max_time = size(rs.seed_log, :timesteps)
decoded_ts = Dict(
    i => ADRIA.analysis.decode_option_ts(
        rs.inputs.option_ts[i], rs.inputs.seed_year_start[i],
        rs.inputs.seed_years[i], rs.inputs.pd_frequency[i], max_time
    )
    for i in 1:length(rs.inputs.option_ts)
)

# ── Filter to parameter subset ────────────────────────────────────────────────

subset_mask =
    rs.inputs.N_seed_CA .== N_seed * N_seed_weights.N_seed_CA .&&
    rs.inputs.N_seed_TA .== N_seed * N_seed_weights.N_seed_TA .&&
    rs.inputs.dhw_scenario .== dhw_scenario .&&
    rs.inputs.min_iv_locations .== min_iv_locations .&&
    rs.inputs.RCP .== parse(Float64, RCP)

subset_idxs = findall(subset_mask)

# ── Find one scenario ID per possible next option ─────────────────────────────

possible_options = [:heat_stress, :geographic_spread, :connectivity,
    :functional_diversity, :balanced]

seed_yr_start = Int(rs.inputs.seed_year_start[subset_idxs[1]])
pd_freq = Int(rs.inputs.pd_frequency[subset_idxs[1]])
n_prev = length(previous_pathway)

function matches_previous(ts)
    return all(
        ts[seed_yr_start + k * pd_freq] == previous_pathway[k + 1]
        for k in 0:(n_prev - 1)
    )
end

selected_scen_ids = Dict{Symbol,Int}()
for next_opt in possible_options
    next_block_start = seed_yr_start + n_prev * pd_freq
    match_idx = findfirst(subset_idxs) do i
        matches_previous(decoded_ts[i]) && decoded_ts[i][next_block_start] == next_opt
    end
    if !isnothing(match_idx)
        selected_scen_ids[next_opt] = subset_idxs[match_idx]
    end
end

# ── Timeframe for the next option block ───────────────────────────────────────

ts_start = seed_yr_start + n_prev * pd_freq   # first timestep of next block
ts_end = ts_start + pd_freq - 1             # last timestep of next block

# ── Compute metrics ───────────────────────────────────────────────────────────

scen_ids = [selected_scen_ids[opt] for opt in possible_options]
n_locs = size(rs.seed_log, :locations)

m_tac = Array(ADRIA.metrics.total_absolute_cover(rs)) .* 1e-6  # km²
loc_hab_area_km2 = rs.loc_area .* rs.loc_max_coral_cover .* 1e-6

# cum_tac: cumulative total absolute cover summed over [ts_start, ts_end]
cum_tac = dropdims(
    sum(m_tac[ts_start:ts_end, :, scen_ids]; dims=1); dims=1
)  # (n_locs, 5)

# cum_rel_tac: cumulative relative cover (TAC / habitat area) over [ts_start, ts_end]
m_rel_tac = m_tac ./ reshape(loc_hab_area_km2, 1, :, 1)
cum_rel_tac = dropdims(
    sum(m_rel_tac[ts_start:ts_end, :, scen_ids]; dims=1); dims=1
)  # (n_locs, 5)

# cum_fd: cumulative coral evenness summed over [ts_start, ts_end]
fd_data = Array(ADRIA.metrics.coral_evenness(rs))  # (timesteps, locs, scenarios)
cum_fd = dropdims(
    sum(fd_data[ts_start:ts_end, :, scen_ids]; dims=1); dims=1
)  # (n_locs, 5)

# ── Comparison: repeat last option ───────────────────────────────────────────
# The comparison is the scenario that continues with previous_pathway[end],
# which is already column cf_col of all metric matrices.

cf_col = findfirst(==(previous_pathway[end]), possible_options)

cum_tac_diff = cum_tac .- cum_tac[:, cf_col]                                          # (n_locs, 5)
cum_rel_tac_diff = cum_rel_tac .- cum_rel_tac[:, cf_col]                                      # (n_locs, 5)
cum_tac_diff_alt = (cum_tac .- cum_tac[:, cf_col]) ./ (cum_tac .+ cum_tac[:, cf_col])         # (n_locs, 5)
cum_fd_diff = cum_fd .- cum_fd[:, cf_col]                                           # (n_locs, 5)
cum_fd_diff_alt = (cum_fd .- cum_fd[:, cf_col]) ./ (cum_fd .+ cum_fd[:, cf_col])          # (n_locs, 5)

# ── Boxplots ──────────────────────────────────────────────────────────────────

scen_labels = string.(possible_options)
x_scens = vcat([fill(i, n_locs) for i in 1:length(possible_options)]...)
cf_label = string(previous_pathway[end])

y_tac = vcat([cum_tac_diff[:, s] for s in 1:length(possible_options)]...)
y_rel_tac = vcat([cum_rel_tac_diff[:, s] for s in 1:length(possible_options)]...)
y_tac_alt = vcat([cum_tac_diff_alt[:, s] for s in 1:length(possible_options)]...)
y_fd = vcat([cum_fd_diff[:, s] for s in 1:length(possible_options)]...)
y_fd_alt = vcat([cum_fd_diff_alt[:, s] for s in 1:length(possible_options)]...)

axis_opts(ylabel) = (;
    xticks=(1:length(possible_options), scen_labels),
    xticklabelrotation=π / 4,
    ylabel
)

fig_box = Figure(; size=(2500, 450))

ax_tac = Axis(fig_box[1, 1]; title="Δ cum. absolute cover vs repeat-$(cf_label)",
    axis_opts("Δ cum. cover [km²]")...)
boxplot!(ax_tac, x_scens, y_tac)

ax_rel_tac = Axis(fig_box[1, 2]; title="Δ cum. relative cover vs repeat-$(cf_label)",
    axis_opts("Δ cum. relative cover")...)
boxplot!(ax_rel_tac, x_scens, y_rel_tac)

ax_tac_alt = Axis(fig_box[1, 3]; title="Normalised Δ cum. cover vs repeat-$(cf_label)",
    axis_opts("(cand−ref)/(cand+ref)")...)
boxplot!(ax_tac_alt, x_scens, y_tac_alt)

ax_fd = Axis(fig_box[1, 4]; title="Δ cum. evenness vs repeat-$(cf_label)",
    axis_opts("Δ cum. evenness")...)
boxplot!(ax_fd, x_scens, y_fd)

ax_fd_alt = Axis(fig_box[1, 5]; title="Normalised Δ cum. evenness vs repeat-$(cf_label)",
    axis_opts("(cand−ref)/(cand+ref)")...)
boxplot!(ax_fd_alt, x_scens, y_fd_alt)

save(
    joinpath(pd_config["plot_output_path"], "benefit_metrics_boxplot.png"),
    fig_box;
    px_per_unit=2
)

# ── Benefit metrics ───────────────────────────────────────────────────────────

"""
Two-sided CVaR benefit mapped to [0,1] via tanh (Metric 1).
When `normalize_by_ref=true` (default), divides the raw CVaR by mean(ref_perf) to make
it dimensionless. Set to false when delta is already dimensionless (e.g. relative cover
or ratio-normalised differences).
Returns 0.5 at break-even, >0.5 net beneficial, <0.5 net detrimental.
"""
function two_sided_cvar(
    delta::AbstractVector, ref_perf::AbstractVector;
    tail_number::Int=150, σ::Float64=0.1, normalize_by_ref::Bool=true
)::Float64
    s = sort(delta)
    cvar_lower = mean(s[1:tail_number])
    cvar_upper = mean(s[(end - tail_number + 1):end])
    raw = cvar_upper - abs(cvar_lower)
    if normalize_by_ref
        ref_mean = mean(ref_perf)
        abs(ref_mean) < eps() && return NaN
        raw = raw / ref_mean
    end
    return (tanh(raw / σ) + 1) / 2
end

"""
Omega ratio (order 1, threshold 0), mapped to [0,1] via tanh(log(omega)) (Metric 2).
Returns 0.5 at break-even, >0.5 net beneficial, <0.5 net detrimental.
"""
function omega_ratio(delta::AbstractVector; σ::Float64=0.1)::Float64
    UPM1 = sum(max(d, 0.0) for d in delta)
    LPM1 = sum(max(-d, 0.0) for d in delta)
    omega = if UPM1 == 0 && LPM1 == 0
        1.0
    elseif LPM1 == 0
        Inf
    else
        UPM1 / LPM1
    end
    log_omega = isinf(omega) ? Inf : log(omega)
    return (tanh(log_omega / σ) + 1) / 2
end

# Parameters to compute benefit
tail_number = 150  # CVaR: number of reefs in each tail
# tanh sensitivity per metric (output ≈ 0.88 when benefit equals σ)
σ_tac = 0.001  # cum_tac_diff       (km², normalized by ref_mean)
σ_rel_tac = 0.01  # cum_rel_tac_diff   (dimensionless relative cover)
σ_tac_alt = 0.01  # cum_tac_diff_alt   ((cand−ref)/(cand+ref))
σ_fd = 0.05  # cum_fd_diff        (evenness units)
σ_fd_alt = 0.001  # cum_fd_diff_alt    ((cand−ref)/(cand+ref))

ref_tac = cum_tac[:, cf_col]
ref_fd = cum_fd[:, cf_col]

for (s, opt) in enumerate(possible_options)
    cvar_tac = two_sided_cvar(
        cum_tac_diff[:, s], ref_tac; tail_number, σ=σ_tac, normalize_by_ref=true
    )
    cvar_rel_tac = two_sided_cvar(
        cum_rel_tac_diff[:, s], ref_tac; tail_number, σ=σ_rel_tac, normalize_by_ref=false
    )
    cvar_tac_alt = two_sided_cvar(
        cum_tac_diff_alt[:, s], ref_tac; tail_number, σ=σ_tac_alt, normalize_by_ref=false
    )
    cvar_fd = two_sided_cvar(
        cum_fd_diff[:, s], ref_fd; tail_number, σ=σ_fd, normalize_by_ref=false
    )
    cvar_fd_alt = two_sided_cvar(
        cum_fd_diff_alt[:, s], ref_fd; tail_number, σ=σ_fd_alt, normalize_by_ref=false
    )

    om_tac = omega_ratio(cum_tac_diff[:, s]; σ=σ_tac)
    om_rel_tac = omega_ratio(cum_rel_tac_diff[:, s]; σ=σ_rel_tac)
    om_tac_alt = omega_ratio(cum_tac_diff_alt[:, s]; σ=σ_tac_alt)
    om_fd = omega_ratio(cum_fd_diff[:, s]; σ=σ_fd)
    om_fd_alt = omega_ratio(cum_fd_diff_alt[:, s]; σ=σ_fd_alt)

    @info "$opt" CVaR_tac = round(cvar_tac; digits=4) CVaR_rel_tac = round(
        cvar_rel_tac; digits=4
    ) CVaR_tac_alt = round(
        cvar_tac_alt; digits=4) CVaR_fd = round(cvar_fd; digits=4) CVaR_fd_alt = round(cvar_fd_alt; digits=4)
end
