#=
Compare the three ResultSets that feed the robustness analyses, to diagnose why they yield
metrics of different magnitude. Loads:

  1. `rs_static`     — a `robustness_analysis_static_options.jl` run. Single ResultSet that
                       already contains its counterfactual as scenario columns (option = -1).
  2. `rs_pathway`    — a pathway simulation as produced by `pd_main.jl`.
  3. `rs_pathway_cf` — the counterfactual as produced by `robustness_analysis_pathways.jl`.

Prints two aligned side-by-side tables: (a) basic data properties (timeframe, dimensions, the
factor values actually simulated, raw cover/evenness magnitudes) and (b) the windowed per-reef
metric summaries — so the three runs line up column-by-column and any mismatch (calendar
timeframe, seed budget, DHW/RCP set, or raw cover magnitude) is obvious at a glance.
=#

include("src/common.jl")

# ── Paths (fill in) ───────────────────────────────────────────────────────────
static_path = ""       # robustness_analysis_static_options.jl run (intervention + counterfactual)
pathway_path = ""      # pd_main.jl pathway simulation
pathway_cf_path = ""   # robustness_analysis_pathways.jl counterfactual

rs_static = ADRIA.load_results(static_path)
rs_pathway = ADRIA.load_results(pathway_path)
rs_pathway_cf = ADRIA.load_results(pathway_cf_path)

# Measurement window (yrs from seed-start), matching both robustness scripts.
horizon = 20

# ── Basic data properties ─────────────────────────────────────────────────────

"""
    rs_facts(result_set, name) -> NamedTuple

Collect timeframe, dimensions, the factor ranges the run actually covers, and raw cover /
evenness magnitudes (whole cube, before any windowing) for one ResultSet.
"""
function rs_facts(result_set, name)
    m_tac = ADRIA.readcubedata(ADRIA.metrics.total_absolute_cover(result_set)) .* 1e-6  # km²
    fd_arr = ADRIA.readcubedata(ADRIA.metrics.coral_evenness(result_set))
    inp = result_set.inputs

    years = try
        collect(lookup(m_tac, :timesteps))
    catch
        nothing
    end
    hab_area_km2 = result_set.loc_area .* result_set.loc_max_coral_cover .* 1e-6
    # GBR-total cover trajectory: sum over locations → (timesteps, scenarios), mean over scens.
    tac_traj = dropdims(sum(m_tac.data; dims=2); dims=2)

    return (
        name=name,
        n_ts=size(m_tac, :timesteps),
        yr_range=isnothing(years) ? "n/a" : "$(first(years))-$(last(years))",
        n_locs=size(m_tac, :locations),
        n_scens=nrow(inp),
        rcps=string(sort(unique(inp.RCP))),
        dhws=string(sort(unique(inp.dhw_scenario))),
        guided=string(sort(unique(inp.guided))),
        min_iv=string(sort(unique(inp.min_iv_locations))),
        nseed_ca=string(round.(extrema(inp.N_seed_CA); sigdigits=3)),
        hab_total=round(sum(hab_area_km2); digits=1),
        cover_mean=round(mean(m_tac.data); digits=4),
        cover_max=round(maximum(m_tac.data); digits=3),
        even_mean=round(mean(fd_arr.data); digits=4),
        gbr_t0=round(mean(tac_traj[1, :]); digits=1),
        gbr_tend=round(mean(tac_traj[end, :]); digits=1)
    )
end

# ── Windowed per-reef metrics (identical formula to both robustness scripts) ──

"""
    reef_metrics(result_set; horizon) -> (n_yrs_above, cum_tac, cum_fd)

Three per-reef metrics over a `horizon`-year window anchored at seed-start, each a
(n_locs, n_scenarios) matrix. Mirrors `robustness_analysis_pathways.jl`.
"""
function reef_metrics(result_set; horizon::Int=horizon)
    n_locs = size(result_set.seed_log, :locations)
    n_scens = nrow(result_set.inputs)

    m_tac = ADRIA.readcubedata(ADRIA.metrics.total_absolute_cover(result_set)) .* 1e-6  # km²
    fd_arr = ADRIA.readcubedata(ADRIA.metrics.coral_evenness(result_set))
    loc_hab_area_km2 = result_set.loc_area .* result_set.loc_max_coral_cover .* 1e-6

    seed_start = Int(result_set.inputs.seed_year_start[1])
    @assert seed_start + horizon - 1 <= size(m_tac, :timesteps) "Run too short for horizon $(horizon)."
    window = seed_start:(seed_start + horizon - 1)

    n_yrs_above = Array{Float64}(undef, n_locs, n_scens)
    for s in 1:n_scens, l in 1:n_locs
        thr = 0.20 * loc_hab_area_km2[l]
        n_yrs_above[l, s] = Float64(
            count(m_tac[timesteps=window, locations=l, scenarios=s].data .>= thr)
        )
    end

    cum_tac = dropdims(sum(m_tac[timesteps=window].data; dims=1); dims=1)
    cum_fd = dropdims(sum(fd_arr[timesteps=window].data; dims=1); dims=1)

    return n_yrs_above, cum_tac, cum_fd
end

"""
    metric_facts(result_set, name; horizon) -> NamedTuple

median / mean / max of each windowed per-reef metric, pooled over reefs × scenarios.
"""
function metric_facts(result_set, name; horizon::Int=horizon)
    nyrs, ctac, cfd = reef_metrics(result_set; horizon=horizon)
    summ(v) = (
        median=round(median(v); digits=3),
        mean=round(mean(v); digits=3),
        max=round(maximum(v); digits=3)
    )
    return (name=name, nyrs=summ(vec(nyrs)), ctac=summ(vec(ctac)), cfd=summ(vec(cfd)))
end

# ── Aligned side-by-side printing ─────────────────────────────────────────────

const _LABELW = 26
const _COLW = 22

_cell(x) = rpad(first(string(x), _COLW - 1), _COLW)

function print_table(title, rows, facts_list)
    println("\n", title)
    println("─"^(_LABELW + _COLW * length(facts_list)))
    print(rpad("", _LABELW))
    for f in facts_list
        print(_cell(f.name))
    end
    println()
    for (label, getter) in rows
        print(rpad(label, _LABELW))
        for f in facts_list
            print(_cell(getter(f)))
        end
        println()
    end
end

# ── Run comparison ────────────────────────────────────────────────────────────

facts = [
    rs_facts(rs_static, "static"),
    rs_facts(rs_pathway, "pathway"),
    rs_facts(rs_pathway_cf, "pathway_cf")
]

property_rows = [
    ("timesteps", f -> f.n_ts),
    ("year range", f -> f.yr_range),
    ("locations", f -> f.n_locs),
    ("scenarios", f -> f.n_scens),
    ("RCPs", f -> f.rcps),
    ("DHW scenarios", f -> f.dhws),
    ("guided", f -> f.guided),
    ("min_iv_locations", f -> f.min_iv),
    ("N_seed_CA range", f -> f.nseed_ca),
    ("habitable km² total", f -> f.hab_total),
    ("cover km²/loc/ts mean", f -> f.cover_mean),
    ("cover km²/loc/ts max", f -> f.cover_max),
    ("evenness mean", f -> f.even_mean),
    ("GBR total cover t0 km²", f -> f.gbr_t0),
    ("GBR total cover tEnd km²", f -> f.gbr_tend)
]

print_table("BASIC DATA PROPERTIES", property_rows, facts)

mfacts = [
    metric_facts(rs_static, "static"; horizon=horizon),
    metric_facts(rs_pathway, "pathway"; horizon=horizon),
    metric_facts(rs_pathway_cf, "pathway_cf"; horizon=horizon)
]

metric_rows = [
    ("n_yrs_above  median", f -> f.nyrs.median),
    ("n_yrs_above  mean", f -> f.nyrs.mean),
    ("n_yrs_above  max", f -> f.nyrs.max),
    ("cum_tac km²·yr median", f -> f.ctac.median),
    ("cum_tac km²·yr mean", f -> f.ctac.mean),
    ("cum_tac km²·yr max", f -> f.ctac.max),
    ("cum_fd  median", f -> f.cfd.median),
    ("cum_fd  mean", f -> f.cfd.mean),
    ("cum_fd  max", f -> f.cfd.max)
]

print_table("WINDOWED PER-REEF METRICS (horizon = $(horizon) yrs)", metric_rows, mfacts)
println()
