#=
Compare the three ResultSets that feed the robustness analyses, to diagnose why they yield
metrics of different magnitude. Loads:

  1. `rs_static`     — a `robustness_analysis_static_options.jl` run. Single ResultSet that
                       already contains its counterfactual as scenario columns (option = -1),
                       unguided (option = 0) and the 5 guided options (option = 1:5).
  2. `rs_pathway`    — a pathway simulation as produced by `pd_main.jl` (`option_ts` per row).
  3. `rs_pathway_cf` — the counterfactual as produced by `robustness_analysis_pathways.jl`.

Analyses
  1. BASIC PROPERTIES + WINDOWED METRICS over every scenario (whole-run overview).
  2. FIXED-PARAMETER subset: restrict all three runs to one (N_seed, min_iv, dhw, rcp) combo
     and compare windowed-metric magnitudes — do they behave similarly under the same params?
  3. OPTION-BY-OPTION: each static guided option vs the pathway whose `option_ts` is constant
     (that option every block). Same params ⇒ expected same/very similar cover & evenness.
  4. COUNTERFACTUAL: static counterfactual (option = -1) vs the pathway counterfactual run.
     Same params ⇒ expected same/very similar cover & evenness.
  5. Δ-DISTRIBUTION: fixed-parameter subset, per-reef Δ against the respective counterfactual
     (pathway − pathway_cf vs static − static_cf) for cumulative cover and cumulative diversity.
     Reports min/max, P10/P90, mean/median of each pooled distribution plus the reef & scenario
     id holding the extreme Δ.

The metrics use one integration window (horizon yrs from seed-start), identical to both
robustness scripts, so magnitudes are apples-to-apples.
=#

include("src/common.jl")

# ── Paths (fill in) ───────────────────────────────────────────────────────────
static_path = ""       # robustness_analysis_static_options.jl run (intervention + counterfactual)
pathway_path = ""      # pd_main.jl pathway simulation
pathway_cf_path = ""   # robustness_analysis_pathways.jl counterfactual

rs_static = ADRIA.load_results(static_path)
rs_pathway = ADRIA.load_results(pathway_path)
rs_pathway_cf = ADRIA.load_results(pathway_cf_path)

# ── Parameters ────────────────────────────────────────────────────────────────
horizon = 20           # measurement window (yrs from seed-start); matches both robustness scripts
pd_frequency = 5       # option-switching period used by the pathway run (pd_main.jl)
N_seed_weights = (
    N_seed_TA=0.15, N_seed_CA=0.5, N_seed_CNA=0.0, N_seed_SM=0.35, N_seed_LM=0.0
)

# Fixed parameter combo for analyses 2–4. Edit to a combo present in all three runs.
sel_N_seed = 1e6
sel_min_iv = 200
sel_dhw = 7
sel_rcp = 45

# ── Load once ─────────────────────────────────────────────────────────────────
# `readcubedata` is the slow step, so each cube is materialised exactly once per ResultSet
# here and reused by every analysis below.

"""
    load_bundle(result_set, name) -> NamedTuple

Materialise the two metric cubes (km²-scaled total absolute cover; coral evenness) once and
bundle them with the ResultSet and per-location habitable area.
"""
function load_bundle(result_set, name)
    m_tac = ADRIA.readcubedata(ADRIA.metrics.total_absolute_cover(result_set)) .* 1e-6  # km²
    fd_arr = ADRIA.readcubedata(ADRIA.metrics.coral_evenness(result_set))
    hab = result_set.loc_area .* result_set.loc_max_coral_cover .* 1e-6  # km²
    return (rs=result_set, name=name, m_tac=m_tac, fd_arr=fd_arr, hab=hab)
end

@info "Loading ResultSets (this materialises the metric cubes; may take a while)…"
static_b = load_bundle(rs_static, "static")
pathway_b = load_bundle(rs_pathway, "pathway")
pathway_cf_b = load_bundle(rs_pathway_cf, "pathway_cf")

option_names = ADRIA.analysis.option_seed_preference().option_name  # 5 guided options, canonical order

# ── Scenario selection helpers ────────────────────────────────────────────────

"""
    param_idxs(rs; N_seed, min_iv, dhw, rcp, match_seed=true) -> Vector{Int}

Row indices of `rs.inputs` matching the given parameters. Each constraint is applied only if
non-`nothing` and the column exists. `match_seed=false` skips the N_seed constraint (for
counterfactual rows that carry no seed budget).
"""
function param_idxs(rs; N_seed=nothing, min_iv=nothing, dhw=nothing, rcp=nothing, match_seed=true)
    inp = rs.inputs
    cols = propertynames(inp)
    mask = trues(nrow(inp))
    if !isnothing(rcp) && (:RCP in cols)
        mask .&= inp.RCP .== float(rcp)
    end
    if !isnothing(dhw) && (:dhw_scenario in cols)
        mask .&= inp.dhw_scenario .== dhw
    end
    if !isnothing(min_iv) && (:min_iv_locations in cols)
        mask .&= inp.min_iv_locations .== min_iv
    end
    if match_seed && !isnothing(N_seed) && (:N_seed_CA in cols)
        mask .&= inp.N_seed_CA .== N_seed * N_seed_weights.N_seed_CA
    end
    return findall(mask)
end

"""
    constant_pathway_option(rs, idx) -> Symbol or nothing

Decode `rs.inputs.option_ts[idx]`; return the single option adopted at every seeding block if
the pathway is constant, otherwise `nothing`.
"""
function constant_pathway_option(rs, idx)
    ss = Int(rs.inputs.seed_year_start[1])
    sy = Int(rs.inputs.seed_years[1])
    mt = size(rs.seed_log, :timesteps)
    ts = ADRIA.analysis.decode_option_ts(rs.inputs.option_ts[idx], ss, sy, pd_frequency, mt; legacy=true)
    window_opts = unique(ts[ss:(ss + sy - 1)])
    return length(window_opts) == 1 ? window_opts[1] : nothing
end

# ── Windowed per-reef metrics (identical window/formula to both robustness scripts) ──

"""
    reef_metrics_sel(bundle, scens; horizon) -> (n_yrs_above, cum_tac, cum_fd)

Three per-reef metrics for the selected scenario indices over a `horizon`-year window anchored
at seed-start, each a (n_locs, n_selected) matrix. Empty selection ⇒ empty matrices.
"""
function reef_metrics_sel(bundle, scens; horizon::Int=horizon)
    isempty(scens) && return (zeros(0, 0), zeros(0, 0), zeros(0, 0))
    seed_start = Int(bundle.rs.inputs.seed_year_start[1])
    nt = size(bundle.m_tac, :timesteps)
    @assert seed_start + horizon - 1 <= nt "$(bundle.name): run too short for horizon $(horizon)."
    window = seed_start:(seed_start + horizon - 1)

    m = bundle.m_tac[timesteps=window, scenarios=scens].data   # (H, locs, nsel)
    fd = bundle.fd_arr[timesteps=window, scenarios=scens].data
    n_locs, nsel = size(m, 2), size(m, 3)

    nyrs = Array{Float64}(undef, n_locs, nsel)
    for s in 1:nsel, l in 1:n_locs
        nyrs[l, s] = count(@view(m[:, l, s]) .>= 0.20 * bundle.hab[l])
    end
    cum_tac = dropdims(sum(m; dims=1); dims=1)   # (locs, nsel)
    cum_fd = dropdims(sum(fd; dims=1); dims=1)
    return nyrs, cum_tac, cum_fd
end

meanor(v) = isempty(v) ? NaN : mean(v)

# ── Aligned printing ──────────────────────────────────────────────────────────

function print_grid(title, colnames, rows; labelw=24, colw=15)
    println("\n", title)
    println("─"^(labelw + colw * length(colnames)))
    print(rpad("", labelw))
    for c in colnames
        print(lpad(first(c, colw - 1), colw))
    end
    println()
    for (label, vals) in rows
        print(rpad(first(label, labelw - 1), labelw))
        for v in vals
            cell = v isa AbstractString ? first(v, colw - 1) : (isnan(v) ? "—" : string(round(v; sigdigits=4)))
            print(lpad(cell, colw))
        end
        println()
    end
end

# ── Analysis 1: basic properties + windowed metrics over all scenarios ────────

"""
    rs_facts(bundle) -> NamedTuple

Timeframe, dimensions, factor coverage and raw (pre-window) cover/evenness magnitudes over
every scenario.
"""
function rs_facts(bundle)
    inp = bundle.rs.inputs
    years = try
        collect(lookup(bundle.m_tac, :timesteps))
    catch
        nothing
    end
    tac_traj = dropdims(sum(bundle.m_tac.data; dims=2); dims=2)  # (timesteps, scenarios)
    fcol(c) = c in propertynames(inp) ? string(sort(unique(inp[!, c]))) : "n/a"
    return (
        name=bundle.name,
        n_ts=size(bundle.m_tac, :timesteps),
        yr_range=isnothing(years) ? "n/a" : "$(first(years))-$(last(years))",
        n_locs=size(bundle.m_tac, :locations),
        n_scens=nrow(inp),
        rcps=fcol(:RCP), dhws=fcol(:dhw_scenario), guided=fcol(:guided), min_iv=fcol(:min_iv_locations),
        nseed_ca=(:N_seed_CA in propertynames(inp)) ? string(round.(extrema(inp.N_seed_CA); sigdigits=3)) : "n/a",
        hab_total=round(sum(bundle.hab); digits=1),
        cover_mean=round(mean(bundle.m_tac.data); digits=4),
        cover_max=round(maximum(bundle.m_tac.data); digits=3),
        even_mean=round(mean(bundle.fd_arr.data); digits=4),
        gbr_t0=round(mean(tac_traj[1, :]); digits=1),
        gbr_tend=round(mean(tac_traj[end, :]); digits=1)
    )
end

facts = [rs_facts(static_b), rs_facts(pathway_b), rs_facts(pathway_cf_b)]
property_rows = [
    ("timesteps", f -> f.n_ts), ("year range", f -> f.yr_range), ("locations", f -> f.n_locs),
    ("scenarios", f -> f.n_scens), ("RCPs", f -> f.rcps), ("DHW scenarios", f -> f.dhws),
    ("guided", f -> f.guided), ("min_iv_locations", f -> f.min_iv), ("N_seed_CA range", f -> f.nseed_ca),
    ("habitable km² total", f -> f.hab_total), ("cover km²/loc/ts mean", f -> f.cover_mean),
    ("cover km²/loc/ts max", f -> f.cover_max), ("evenness mean", f -> f.even_mean),
    ("GBR cover t0 km²", f -> f.gbr_t0), ("GBR cover tEnd km²", f -> f.gbr_tend)
]
print_grid("[1] BASIC DATA PROPERTIES (all scenarios)",
    [f.name for f in facts],
    [(lab, [string(get(f)) for f in facts]) for (lab, get) in property_rows]; labelw=24, colw=18)

"""
    metric_summary(bundle, scens) -> NamedTuple of (median, mean, max) per metric.
"""
function metric_summary(bundle, scens)
    nyrs, ctac, cfd = reef_metrics_sel(bundle, scens)
    s(v) = (median=median(vec(v)), mean=mean(vec(v)), max=maximum(vec(v)))
    isempty(scens) && return (nyrs=(median=NaN, mean=NaN, max=NaN), ctac=(median=NaN, mean=NaN, max=NaN), cfd=(median=NaN, mean=NaN, max=NaN))
    return (nyrs=s(nyrs), ctac=s(ctac), cfd=s(cfd))
end

all_scens(b) = 1:nrow(b.rs.inputs)
msum = [metric_summary(static_b, all_scens(static_b)),
    metric_summary(pathway_b, all_scens(pathway_b)),
    metric_summary(pathway_cf_b, all_scens(pathway_cf_b))]
metric_rows = [
    ("n_yrs_above  median", m -> m.nyrs.median), ("n_yrs_above  mean", m -> m.nyrs.mean), ("n_yrs_above  max", m -> m.nyrs.max),
    ("cum_tac km²·yr median", m -> m.ctac.median), ("cum_tac km²·yr mean", m -> m.ctac.mean), ("cum_tac km²·yr max", m -> m.ctac.max),
    ("cum_fd  median", m -> m.cfd.median), ("cum_fd  mean", m -> m.cfd.mean), ("cum_fd  max", m -> m.cfd.max)
]
print_grid("[1] WINDOWED PER-REEF METRICS (horizon = $(horizon) yrs, all scenarios)",
    ["static", "pathway", "pathway_cf"],
    [(lab, [get(m) for m in msum]) for (lab, get) in metric_rows])

# ── Analysis 2: fixed-parameter subset ────────────────────────────────────────

sel_str = "N_seed=$(sel_N_seed), min_iv=$(sel_min_iv), dhw=$(sel_dhw), rcp=$(sel_rcp)"

stat_sel = param_idxs(static_b.rs; N_seed=sel_N_seed, min_iv=sel_min_iv, dhw=sel_dhw, rcp=sel_rcp)
path_sel = param_idxs(pathway_b.rs; N_seed=sel_N_seed, min_iv=sel_min_iv, dhw=sel_dhw, rcp=sel_rcp)
path_cf_sel = param_idxs(pathway_cf_b.rs; dhw=sel_dhw, rcp=sel_rcp, match_seed=false)

@info "[2] matched scenarios" static=length(stat_sel) pathway=length(path_sel) pathway_cf=length(path_cf_sel)

msum2 = [metric_summary(static_b, stat_sel),
    metric_summary(pathway_b, path_sel),
    metric_summary(pathway_cf_b, path_cf_sel)]
print_grid("[2] WINDOWED METRICS — fixed params ($(sel_str))",
    ["static", "pathway", "pathway_cf"],
    [(lab, [get(m) for m in msum2]) for (lab, get) in metric_rows])

# ── Analysis 3: static option vs constant pathway, same params ────────────────

# Pathway scenarios matching the fixed params, grouped by their constant starting/only option.
path_const = Dict{Symbol,Vector{Int}}()
for i in path_sel
    o = constant_pathway_option(pathway_b.rs, i)
    isnothing(o) && continue
    push!(get!(path_const, o, Int[]), i)
end

opt3_rows = Vector{Tuple{String,Vector{Any}}}()
for option in option_names
    o = findfirst(==(option), option_names)
    stat_o = intersect(stat_sel, findall(static_b.rs.inputs.option .== o))
    path_o = get(path_const, option, Int[])

    _, stac, sfd = reef_metrics_sel(static_b, stat_o)
    _, ptac, pfd = reef_metrics_sel(pathway_b, path_o)
    st, pt = meanor(vec(stac)), meanor(vec(ptac))
    sf, pf = meanor(vec(sfd)), meanor(vec(pfd))
    push!(opt3_rows, (string(option), Any[st, pt, pt / st, sf, pf, pf / sf]))
end
@info "[3] pathway constant-option scenario counts" (; (o => length(get(path_const, o, Int[])) for o in option_names)...)
print_grid("[3] STATIC OPTION vs CONSTANT PATHWAY — mean over reefs ($(sel_str))",
    ["cum_tac stat", "cum_tac path", "tac ratio", "cum_fd stat", "cum_fd path", "fd ratio"],
    opt3_rows; labelw=22, colw=15)

# ── Analysis 4: static counterfactual vs pathway counterfactual, same params ──

stat_cf = intersect(param_idxs(static_b.rs; min_iv=sel_min_iv, dhw=sel_dhw, rcp=sel_rcp, match_seed=false),
    findall(static_b.rs.inputs.option .== -1))
# pathway_cf run is entirely counterfactual; match on params it carries.
path_cf = path_cf_sel

@info "[4] counterfactual scenario counts" static=length(stat_cf) pathway_cf=length(path_cf)

_, scf_tac, scf_fd = reef_metrics_sel(static_b, stat_cf)
_, pcf_tac, pcf_fd = reef_metrics_sel(pathway_cf_b, path_cf)
st_tac, pt_tac = meanor(vec(scf_tac)), meanor(vec(pcf_tac))
st_fd, pt_fd = meanor(vec(scf_fd)), meanor(vec(pcf_fd))
print_grid("[4] COUNTERFACTUAL: static (option=-1) vs pathway_cf — mean over reefs ($(sel_str))",
    ["static", "pathway_cf", "ratio"],
    [("cum_tac km²·yr", Any[st_tac, pt_tac, pt_tac / st_tac]),
        ("cum_fd", Any[st_fd, pt_fd, pt_fd / st_fd])]; labelw=18, colw=16)
println()

# ── Analysis 5: distribution of Δ-vs-counterfactual, pathway vs static ────────
#
# Same fixed-parameter subset and horizon window as analysis 2. For each performance metric
# (cumulative cover cum_tac; cumulative diversity cum_fd) build the per-reef Δ against the
# *respective* counterfactual:
#     pathway Δ = pathway_metric[reef, scen]  − mean(pathway_cf) over reef
#     static  Δ = static_metric [reef, scen]  − mean(static_cf)  over reef
# then compare the two pooled Δ distributions (over reefs × matched scenarios): min/max,
# P10/P90, mean/median, plus the (reef, scenario id) carrying the extreme Δ.

"""
    loc_labels(bundle) -> Vector

Location identifiers along the metric cube's `:locations` dim (falls back to 1-based indices).
"""
function loc_labels(bundle)
    try
        return collect(lookup(bundle.m_tac, :locations))
    catch
        return collect(1:size(bundle.m_tac, :locations))
    end
end

"""
    cf_vector(bundle, cf_scens) -> (ctac, cfd)

Per-reef counterfactual cum_tac / cum_fd from the first matched counterfactual scenario
(matches the `findfirst` convention in the robustness scripts).
"""
function cf_vector(bundle, cf_scens)
    _, ctac, cfd = reef_metrics_sel(bundle, cf_scens[1:1])
    return (ctac=vec(ctac), cfd=vec(cfd))
end

"""
    delta_mats(bundle, scens, cf_vec) -> (ctac, cfd)

Per-reef Δ (intervention − counterfactual) matrices, each (n_locs, n_scens). The cf vectors
broadcast column-wise (subtracted from every scenario).
"""
function delta_mats(bundle, scens, cf_vec)
    _, ctac, cfd = reef_metrics_sel(bundle, scens)
    return (ctac=ctac .- cf_vec.ctac, cfd=cfd .- cf_vec.cfd)
end

"""
    dist_stats(delta) -> NamedTuple

Pooled distribution summary (min, P10, median, mean, P90, max) over reefs × scenarios.
"""
function dist_stats(delta)
    v = vec(delta)
    return (
        min=minimum(v), p10=quantile(v, 0.10), median=median(v),
        mean=mean(v), p90=quantile(v, 0.90), max=maximum(v)
    )
end

"""
    extreme_info(delta, scens, labels) -> NamedTuple

(reef, scenario id, value) for the min and max Δ in a (n_locs, n_scens) matrix. `scens` maps
local scenario columns back to global input rows; `labels` maps local rows to location ids.
"""
function extreme_info(delta, scens, labels)
    imin, imax = argmin(delta), argmax(delta)
    return (
        min_reef=labels[imin[1]], min_scen=scens[imin[2]], min_val=delta[imin],
        max_reef=labels[imax[1]], max_scen=scens[imax[2]], max_val=delta[imax]
    )
end

pathway_cf_vec = cf_vector(pathway_cf_b, path_cf_sel)
static_cf_vec = cf_vector(static_b, stat_cf)

path_delta = delta_mats(pathway_b, path_sel, pathway_cf_vec)
stat_delta = delta_mats(static_b, stat_sel, static_cf_vec)

path_labels = loc_labels(pathway_b)
stat_labels = loc_labels(static_b)

# Columns: pathway/static × (cover, diversity). Rows: distribution summary.
d5 = (
    path_tac=dist_stats(path_delta.ctac), stat_tac=dist_stats(stat_delta.ctac),
    path_fd=dist_stats(path_delta.cfd), stat_fd=dist_stats(stat_delta.cfd)
)
stat_rows = [
    ("min", s -> s.min), ("P10", s -> s.p10), ("median", s -> s.median),
    ("mean", s -> s.mean), ("P90", s -> s.p90), ("max", s -> s.max)
]
print_grid("[5] Δ vs COUNTERFACTUAL distribution — $(sel_str)",
    ["cover path", "cover static", "div path", "div static"],
    [(lab, [get(d5.path_tac), get(d5.stat_tac), get(d5.path_fd), get(d5.stat_fd)])
     for (lab, get) in stat_rows]; labelw=12, colw=15)

# Reef + scenario id carrying the extreme Δ for each distribution.
extremes = [
    ("cover  pathway", extreme_info(path_delta.ctac, path_sel, path_labels)),
    ("cover  static", extreme_info(stat_delta.ctac, stat_sel, stat_labels)),
    ("div    pathway", extreme_info(path_delta.cfd, path_sel, path_labels)),
    ("div    static", extreme_info(stat_delta.cfd, stat_sel, stat_labels))
]
println("\n[5] EXTREME Δ locations & scenario ids — $(sel_str)")
println("─"^96)
println(rpad("distribution", 16), rpad("min Δ (reef @ scenario)", 40), "max Δ (reef @ scenario)")
for (name, e) in extremes
    minstr = "$(round(e.min_val; sigdigits=4))  ($(e.min_reef) @ scen $(e.min_scen))"
    maxstr = "$(round(e.max_val; sigdigits=4))  ($(e.max_reef) @ scen $(e.max_scen))"
    println(rpad(name, 16), rpad(minstr, 40), maxstr)
end
println()

# ── Analysis 6: confirm the pathway Δ is dominated by unpaired cross-run noise ─
#
# Test 1 ruled out location misordering (identical_order = true) and Test 2 showed two
# *independent* no-intervention runs already differ per reef (P10 −0.01, P90 +0.089) — the
# same order as the pathway "intervention" Δ. The block below confirms hypothesis #2: the
# static Δ is tight only because it is a *paired* within-ResultSet comparison, while the
# pathway Δ subtracts an independently-sampled counterfactual whose run-to-run differences
# swamp the intervention signal.

# Align a per-reef vector expressed in ResultSet `from_b`'s location order onto `to_b`'s order,
# matching by location id (Test 1 says the sets agree; this makes the code order-independent).
function align_by_id(vals, from_b, to_b)
    pos = Dict(id => i for (i, id) in enumerate(loc_labels(from_b)))
    to_ids = loc_labels(to_b)
    @assert issubset(Set(to_ids), Set(keys(pos))) "location id sets differ between ResultSets"
    return [vals[pos[id]] for id in to_ids]
end

# ── (6a) Factor-pairing audit — which non-intervention factors vary within a selection ──
# If none vary across the static intervention+cf rows, those rows are fully paired (identical
# environment/free factors) and the tiny static Δ is the genuine intervention effect.
iv_cols = Set([:option, :guided, :N_seed_TA, :N_seed_CA, :N_seed_CNA, :N_seed_SM, :N_seed_LM, :option_ts])

function varying_factors(inp, rows)
    out = Pair{Symbol,Any}[]
    for c in propertynames(inp)
        c in iv_cols && continue
        vals = unique(inp[rows, c])
        length(vals) > 1 && push!(out, c => vals)
    end
    return out
end

@info "[6a] static — non-intervention factors varying across intervention+cf rows (empty ⇒ fully paired)" varying_factors(static_b.rs.inputs, vcat(stat_sel, stat_cf))
@info "[6a] pathway — non-intervention factors varying across intervention rows" varying_factors(pathway_b.rs.inputs, path_sel)

# Does the single pathway_cf row share each factor's value with the pathway intervention rows?
cf_row = path_cf_sel[1]
cf_mismatches = Symbol[]
for c in propertynames(pathway_b.rs.inputs)
    (c in iv_cols) && continue
    (c in propertynames(pathway_cf_b.rs.inputs)) || continue
    pvals = unique(pathway_b.rs.inputs[path_sel, c])
    (length(pvals) == 1 && pathway_cf_b.rs.inputs[cf_row, c] == pvals[1]) || push!(cf_mismatches, c)
end
@info "[6a] pathway_cf factors NOT matching the intervention (each ⇒ unpaired noise source)" cf_mismatches

# ── (6b) Noise floor vs signal — three Δ distributions, same metric, side by side ──
# cf-vs-cf is the pure unpaired no-intervention baseline (Test 2), aligned onto static reefs.
pathway_cf_on_static = (
    ctac=align_by_id(pathway_cf_vec.ctac, pathway_cf_b, static_b),
    cfd=align_by_id(pathway_cf_vec.cfd, pathway_cf_b, static_b)
)
cfcf = (
    ctac=static_cf_vec.ctac .- pathway_cf_on_static.ctac,
    cfd=static_cf_vec.cfd .- pathway_cf_on_static.cfd
)
print_grid("[6b] Δ DISTRIBUTIONS — pathway signal vs static signal vs unpaired noise floor",
    ["cov path Δ", "cov stat Δ", "cov cf−cf", "div path Δ", "div stat Δ", "div cf−cf"],
    [(lab, [get(dist_stats(path_delta.ctac)), get(dist_stats(stat_delta.ctac)), get(dist_stats(cfcf.ctac)),
        get(dist_stats(path_delta.cfd)), get(dist_stats(stat_delta.cfd)), get(dist_stats(cfcf.cfd))])
     for (lab, get) in stat_rows]; labelw=10, colw=14)

# ── (6c) Controlled experiment — same intervention numerator, swap paired→unpaired cf ──
# Static intervention minus its OWN paired cf (tight) vs minus the independent pathway_cf
# (should inflate to pathway-like spread). If it does, the pathway Δ width is the pairing
# artefact, not a stronger intervention.
_, stat_int_ctac, stat_int_cfd = reef_metrics_sel(static_b, stat_sel)
stat_unpaired = (
    ctac=stat_int_ctac .- pathway_cf_on_static.ctac,
    cfd=stat_int_cfd .- pathway_cf_on_static.cfd
)
print_grid("[6c] STATIC intervention Δ — paired cf vs independent (unpaired) cf",
    ["cov paired", "cov unpair", "div paired", "div unpair"],
    [(lab, [get(dist_stats(stat_delta.ctac)), get(dist_stats(stat_unpaired.ctac)),
        get(dist_stats(stat_delta.cfd)), get(dist_stats(stat_unpaired.cfd))])
     for (lab, get) in stat_rows]; labelw=10, colw=14)
println()

# ── Analysis 7: is MCB enabled in each ResultSet, and in how many scenarios ────
#
# MCB modifies the DHW environment only when BOTH mcb_albedo > 0 and mcb_duration > 0
# (scenario.jl gate `mcb_duration > 0.0 && mcb_albedo > 0.0`); dhw_scenario is always > 0 in
# these runs and mcb_deployment_freq only sets the cadence, so those two factors alone decide
# whether MCB happens. A ResultSet where MCB is active is not comparable, per reef, to one
# where it is off — a candidate source of the cross-run Δ noise seen above.

"""
    mcb_status(bundle) -> NamedTuple

Whether MCB is switchable in this ResultSet (both factors present) and, if so, how many
scenarios have it active plus the value ranges. `has_cols=false` ⇒ MCB machinery inactive.
"""
function mcb_status(bundle)
    inp = bundle.rs.inputs
    cols = propertynames(inp)
    (:mcb_albedo in cols && :mcb_duration in cols) || return (name=bundle.name, has_cols=false)
    albedo, duration = inp.mcb_albedo, inp.mcb_duration
    active = (albedo .> 0.0) .& (duration .> 0.0)
    return (
        name=bundle.name, has_cols=true, n=nrow(inp), n_active=count(active),
        frac_active=mean(active), albedo_rng=extrema(albedo), duration_rng=extrema(duration),
        albedo_vals=sort(unique(albedo)), duration_vals=sort(unique(duration))
    )
end

println("\n[7] MCB ACTIVITY per ResultSet (active ⇔ mcb_albedo>0 AND mcb_duration>0)")
println("─"^96)
println(rpad("resultset", 12), rpad("scenarios", 11), rpad("MCB active", 20),
    rpad("albedo range", 26), "duration range")
for b in (static_b, pathway_b, pathway_cf_b)
    s = mcb_status(b)
    if !s.has_cols
        println(rpad(s.name, 12), "MCB factors absent from inputs — MCB machinery inactive")
        continue
    end
    println(
        rpad(s.name, 12),
        rpad(string(s.n), 11),
        rpad("$(s.n_active) ($(round(100 * s.frac_active; digits=1))%)", 20),
        rpad(string(s.albedo_rng), 26),
        string(s.duration_rng)
    )
end
for b in (static_b, pathway_b, pathway_cf_b)
    s = mcb_status(b)
    s.has_cols && @info "[7] $(s.name) MCB factor values" albedo=s.albedo_vals duration=s.duration_vals
end
println()
