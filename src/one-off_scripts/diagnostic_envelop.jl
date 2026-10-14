#=
Diagnostics for the outcome-envelope figure (`plot_envelope_timeseries`).

Puzzle this file investigates
-----------------------------
In the envelope figure the *starting option* (decision block 1) clearly separates the median
curves, yet the 125 downstream pathways *within* a starting option produce an almost-collapsed
P10–P90 band — and the collapse is worst for cumulative cover at the low seed budget (1e6), while
the same setting at 1e8 shows a visible band. We want to know whether the downstream option
choices (blocks 2–4) genuinely do little, or whether the analysis hides a real effect.

Domain facts (from the modeller)
--------------------------------
- 3806 reefs total; tail_number (150) selects the best/worst tail.
- The seed *count* does not feed reef selection directly, but seeding changes the decision-matrix
  values (cover, depth, connectivity, evenness over time), so the budget can shift selection
  *indirectly*. Hence the seeded-reef set (and the tail) is NOT assumed budget-invariant — it is
  measured here.

Hypotheses
----------
  H1  Starting-option-locked reef set: within a starting option all 125 pathways seed ~the same
      reefs (same tail); different starting options seed different reefs.
  H2  Early-cohort / magnitude dominance: block-1 corals grow the whole window and dominate cover;
      blocks 2–4 add little in absolute terms at 1e6, visibly more at 1e8. (Downstream effect is
      real but sub-threshold at low budget.)
  H3  Shared-refugia, magnitude-only: the tail is ~the same reefs even across starting options;
      separation is only in how much cover those same reefs gain.
  H4  Deployment insensitive to downstream blocks: seed_log locations barely change with the
      downstream option, so no metric could show within-option spread.
  H5  Metric shape: cover (a compounding stock) is downstream-insensitive; evenness / years-above
      retain within-option spread.

Assumptions about scope (IMPORTANT)
-----------------------------------
This script is meant to be run *after* `src/robustness_analysis.jl`, so every argument of
`envelope_timeseries` is already defined in `Main`, namely:

    rs, option_names, opt_metric_mats, base,
    dhw_scenarios, tail_number, seed_years, horizon,
    seed_year_start, pd_frequency, N_seed_weights, ps

plus the data-processing helpers it calls (all from `src/data_processing/robustness.jl`, already
`include`d by the main script): `pathway_dhw_scenarios`, `delta_tail_ratio`, `tail_indices`,
`find_counterfactual`, `tail_series_vs_cf`, `compact_pathway`.

The two budgets compared below are taken from `ps.n_locations` with N_seed ∈ {1e6, 1e8}; edit
`DIAG_BUDGETS` / `DIAG_NLOC` if your run uses different values.
=#

using Statistics, StatsBase

# ── Guard: fail loudly if the expected context is missing ─────────────────────
for nm in (:rs, :option_names, :opt_metric_mats, :base, :dhw_scenarios, :tail_number,
    :seed_years, :horizon, :seed_year_start, :pd_frequency, :N_seed_weights, :ps)
    @assert isdefined(Main, nm) "`$nm` not in scope — run src/robustness_analysis.jl first."
end

# ── Diagnostic configuration ──────────────────────────────────────────────────
DIAG_BUDGETS = [1e6, 1e8]        # compare low vs high seed budget
DIAG_NLOC = ps.n_locations       # keep the n_locations of the figure under investigation
DIAG_DHW = (7 in dhw_scenarios) ? 7 : first(dhw_scenarios)  # fixed DHW for block decomposition
JACC_HIGH = 0.90                 # "sets essentially coincide" threshold
ETA2_SMALL = 0.05                # "block explains ~no variance" threshold

# Focus metric for most tests: cumulative cover (idx 2 in opt_metric_mats, sym :tac), best-150 tail.
const METRICS = (nyrs=(idx=1, sym=:nyrs), tac=(idx=2, sym=:tac), fd=(idx=3, sym=:fd))

_norm(sym) = sym === :nyrs ? float(seed_years) : 0.0
jaccard(a, b) = (u = length(union(a, b)); u == 0 ? 1.0 : length(intersect(a, b)) / u)

"""Per-location set (Vector{Int}) of reefs ever seeded in scenario column `s` (any timestep/coral)."""
function seeded_locations(s::Int)
    tot = sum(rs.seed_log[scenarios=s]; dims=(:timesteps, :coral_id))
    return findall(vec(tot.data) .> 0)
end

"""
    select_pathways(N_seed; metric, side, col) -> Dict(option => Dict(option_ts => rec))

For parameter set `(N_seed, DIAG_NLOC)`, replicate `envelope_timeseries`' per-pathway selection for
one metric/side/column: choose each pathway's best/worst DHW, its tail, and the plotted curve.
`rec` = (option_ts, dhw, s, tail, val, final, curve).
"""
function select_pathways(N_seed; metric=METRICS.tac, side::Symbol=:best, col::Symbol=:best)
    scalar_mat = opt_metric_mats[metric.idx]
    m_tac, fd_arr, hab, seed_start = base
    window = seed_start:(seed_start + horizon - 1)
    norm = _norm(metric.sym)
    cf_col = Dict(d => find_counterfactual(rs, d) for d in dhw_scenarios)

    paths = pathway_dhw_scenarios(
        rs, option_names;
        ps=(N_seed=N_seed, n_locations=DIAG_NLOC), dhw_scenarios,
        seed_year_start, seed_years, pd_frequency, N_seed_weights
    )

    out = Dict{eltype(option_names),Dict{Any,NamedTuple}}()
    for option in option_names
        recs = Dict{Any,NamedTuple}()
        for (key, dhw_map) in paths[option]
            best = nothing
            for (dhw, s) in dhw_map
                b, t, _ = delta_tail_ratio(scalar_mat[:, s], scalar_mat[:, cf_col[dhw]];
                    tail_number, norm)
                val = side === :worst ? b : t
                if isnothing(best) || (col === :worst ? val < best.val : val > best.val)
                    bi, ti = tail_indices(scalar_mat[:, s] .- scalar_mat[:, cf_col[dhw]]; tail_number)
                    best = (; val, dhw, s, tail=(side === :worst ? bi : ti))
                end
            end
            curve = tail_series_vs_cf(m_tac, fd_arr, hab, window,
                best.s, cf_col[best.dhw], best.tail, metric.sym)
            recs[key] = (; option_ts=key, best.dhw, best.s, best.tail, best.val,
                final=curve[end], curve)
        end
        out[option] = recs
    end
    return out
end

# Mean pairwise Jaccard of a list of index-sets (all pairs; 125C2 is cheap).
function mean_pairwise_jaccard(sets)
    n = length(sets)
    n < 2 && return NaN
    acc = 0.0; cnt = 0
    for i in 1:(n - 1), j in (i + 1):n
        acc += jaccard(sets[i], sets[j]); cnt += 1
    end
    return acc / cnt
end

# Representative (constant-sequence) pathway of a starting option: the one that adopts `option`
# in every block. `recs` is the inner Dict from `select_pathways`.
function representative_rec(option, recs)
    mt = size(rs.seed_log, :timesteps)
    for (key, rec) in recs
        seq = compact_pathway(key, seed_year_start, seed_years, pd_frequency, mt)
        all(==(option), seq) && return rec
    end
    return first(values(recs))  # fallback
end

# η² (fraction of variance explained) of grouping `x` by `labels`.
function eta_squared(x::AbstractVector, labels::AbstractVector)
    grand = mean(x)
    ss_tot = sum((x .- grand) .^ 2)
    ss_tot == 0 && return 0.0
    ss_between = 0.0
    for g in unique(labels)
        xs = x[labels .== g]
        ss_between += length(xs) * (mean(xs) - grand)^2
    end
    return ss_between / ss_tot
end

# ── TEST 1 — tail-set Jaccard: within option, between options, across budgets ──
function test_tail_jaccard()
    println("\n=== TEST 1: tail-set overlap (cumulative cover, best-150) ===")
    sel = Dict(N => select_pathways(N) for N in DIAG_BUDGETS)

    within = Dict{Float64,Dict}()
    for N in DIAG_BUDGETS
        within[N] = Dict(o => mean_pairwise_jaccard([r.tail for r in values(sel[N][o])])
                         for o in option_names)
        println("  N=$(N): within-option mean pairwise tail Jaccard")
        for o in option_names
            println("      $(rpad(string(o), 22)) $(round(within[N][o]; digits=3))  " *
                    "(n=$(length(sel[N][o])))")
        end
        reps = [representative_rec(o, sel[N][o]).tail for o in option_names]
        println("      between-option (pure pathways): $(round(mean_pairwise_jaccard(reps); digits=3))")
    end

    # Across budgets: same option_ts, compare tails.
    println("  across-budget tail Jaccard (same option_ts, per option):")
    lo, hi = DIAG_BUDGETS[1], DIAG_BUDGETS[end]
    for o in option_names
        keys_both = intersect(keys(sel[lo][o]), keys(sel[hi][o]))
        js = [jaccard(sel[lo][o][k].tail, sel[hi][o][k].tail) for k in keys_both]
        println("      $(rpad(string(o), 22)) $(round(mean(js); digits=3))")
    end
    return sel, within
end

# ── TEST 2 — seed_log location sensitivity to downstream blocks ────────────────
function test_seedlog_sensitivity()
    println("\n=== TEST 2: seeded-location set vs downstream option (DHW=$(DIAG_DHW)) ===")
    paths = pathway_dhw_scenarios(
        rs, option_names;
        ps=(N_seed=DIAG_BUDGETS[1], n_locations=DIAG_NLOC), dhw_scenarios,
        seed_year_start, seed_years, pd_frequency, N_seed_weights
    )
    reps = NamedTuple[]
    for o in option_names
        sets = Vector{Int}[]
        for (key, dhw_map) in paths[o]
            haskey(dhw_map, DIAG_DHW) || continue
            push!(sets, seeded_locations(dhw_map[DIAG_DHW]))
        end
        wj = mean_pairwise_jaccard(sets)
        println("  $(rpad(string(o), 22)) within-option seeded-set Jaccard = $(round(wj; digits=3))  " *
                "(n=$(length(sets)), |set|≈$(isempty(sets) ? 0 : round(Int, mean(length.(sets)))))")
        push!(reps, (; o, set=isempty(sets) ? Int[] : sets[1]))
    end
    bj = mean_pairwise_jaccard([r.set for r in reps])
    println("  between starting options (one pathway each): $(round(bj; digits=3))")
    return nothing
end

# ── TEST 3 — variance of the plotted endpoint explained by each block ──────────
function test_block_variance(; N=DIAG_BUDGETS[1], metric=METRICS.tac, side=:best)
    println("\n=== TEST 3: endpoint-Δ variance by decision block " *
            "(N=$(N), metric=$(metric.sym), side=$(side), DHW=$(DIAG_DHW)) ===")
    m_tac, fd_arr, hab, seed_start = base
    window = seed_start:(seed_start + horizon - 1)
    mt = size(rs.seed_log, :timesteps)
    scalar_mat = opt_metric_mats[metric.idx]

    paths = pathway_dhw_scenarios(
        rs, option_names;
        ps=(N_seed=N, n_locations=DIAG_NLOC), dhw_scenarios,
        seed_year_start, seed_years, pd_frequency, N_seed_weights
    )
    cf = find_counterfactual(rs, DIAG_DHW)

    option_medians = Float64[]
    within_var = Float64[]
    for o in option_names
        finals = Float64[]; blk = Dict(b => Symbol[] for b in 2:4)
        for (key, dhw_map) in paths[o]
            haskey(dhw_map, DIAG_DHW) || continue
            s = dhw_map[DIAG_DHW]
            bi, ti = tail_indices(scalar_mat[:, s] .- scalar_mat[:, cf]; tail_number)
            tail = side === :worst ? bi : ti
            curve = tail_series_vs_cf(m_tac, fd_arr, hab, window, s, cf, tail, metric.sym)
            push!(finals, curve[end])
            seq = compact_pathway(key, seed_year_start, seed_years, pd_frequency, mt)
            for b in 2:4
                push!(blk[b], seq[b])
            end
        end
        isempty(finals) && continue
        push!(option_medians, median(finals))
        push!(within_var, var(finals))
        cv = std(finals) / abs(median(finals) + eps())
        etas = join(["blk$(b)=$(round(eta_squared(finals, blk[b]); digits=3))" for b in 2:4], "  ")
        println("  $(rpad(string(o), 22)) median=$(round(median(finals); sigdigits=3)) " *
                "CV=$(round(cv; digits=3))  η²: $etas")
    end
    bv = var(option_medians)           # between-starting-option variance of medians
    wv = mean(within_var)              # mean within-option variance
    println("  between-starting-option var(median) = $(round(bv; sigdigits=3))")
    println("  mean within-option var              = $(round(wv; sigdigits=3))")
    println("  within / between ratio              = $(round(wv / (bv + eps()); sigdigits=3))")
    return (; bv, wv)
end

# ── TEST 4 — budget scaling of within-option spread (CV) ───────────────────────
function test_budget_scaling()
    println("\n=== TEST 4: within-option CV of endpoint-Δ, 1e6 vs 1e8 (cover, best-150) ===")
    for N in DIAG_BUDGETS
        sel = select_pathways(N)
        cvs = Float64[]
        for o in option_names
            finals = [r.final for r in values(sel[o])]
            length(finals) < 2 && continue
            push!(cvs, std(finals) / abs(median(finals) + eps()))
        end
        println("  N=$(N): median across options of within-option CV = $(round(median(cvs); digits=3))")
    end
    return nothing
end

# ── TEST 5 — metric contrast at 1e6 (does cover collapse while others don't?) ──
function test_metric_contrast()
    println("\n=== TEST 5: block-2 η² at 1e6 by metric (downstream signal per metric) ===")
    for (name, metric) in pairs(METRICS)
        res = test_block_variance(; N=DIAG_BUDGETS[1], metric, side=:best)  # prints its own lines
        println("    → metric=$(name): within/between = $(round(res.wv / (res.bv + eps()); sigdigits=3))")
    end
    return nothing
end

# ── Run + verdict ─────────────────────────────────────────────────────────────
sel1, within1 = test_tail_jaccard()
test_seedlog_sensitivity()
cover_var = test_block_variance()
test_budget_scaling()

println("\n========================= LIKELY EXPLANATION =========================")
# Summaries used for the automated verdict (cover, best-150).
within_tail_lo = median(collect(values(within1[DIAG_BUDGETS[1]])))
# Re-measure seeded-set within-option Jaccard (cover-independent) for the verdict.
let
    paths = pathway_dhw_scenarios(
        rs, option_names;
        ps=(N_seed=DIAG_BUDGETS[1], n_locations=DIAG_NLOC), dhw_scenarios,
        seed_year_start, seed_years, pd_frequency, N_seed_weights
    )
    seed_js = Float64[]
    for o in option_names
        sets = [seeded_locations(dm[DIAG_DHW]) for (_, dm) in paths[o] if haskey(dm, DIAG_DHW)]
        push!(seed_js, mean_pairwise_jaccard(sets))
    end
    seed_within = median(seed_js)

    println("within-option tail Jaccard (cover) ≈ $(round(within_tail_lo; digits=3))")
    println("within-option seeded-set Jaccard   ≈ $(round(seed_within; digits=3))")
    println("within/between variance ratio      ≈ $(round(cover_var.wv / (cover_var.bv + eps()); sigdigits=3))")
    println()

    if seed_within > JACC_HIGH && within_tail_lo > JACC_HIGH
        println("→ H1/H4 favoured: downstream blocks barely change WHERE seeding happens, so the")
        println("  starting option effectively locks the reef set and the band must collapse for")
        println("  every metric. Check TEST 5 — if evenness/years also collapse, H4; else H1.")
    elseif seed_within <= JACC_HIGH
        println("→ H2/H5 favoured: downstream blocks DO change the seeded reefs (seeded-set Jaccard")
        println("  below $(JACC_HIGH)), but the cover consequence is tiny at 1e6. If TEST 4 shows CV")
        println("  rising with budget, and TEST 5 shows cover η² ≈ 0 while evenness/years η² > 0,")
        println("  the collapse is a magnitude / metric-shape effect (H2 + H5), not hidden-by-analysis.")
        if within_tail_lo > JACC_HIGH
            println("  (Tail Jaccard still high despite differing seeded sets ⇒ the tail is a stable")
            println("   near-refugia set: H3 — separation is magnitude-only on shared reefs.)")
        end
    else
        println("→ Mixed signal: inspect the per-test tables above; no single hypothesis dominates.")
    end
end
println("Run test_metric_contrast() for the per-metric η² breakdown (TEST 5).")
println("======================================================================")
