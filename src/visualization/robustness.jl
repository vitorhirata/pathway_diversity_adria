#=
Visualisation for the robustness analyses.

Plotting functions for both `robustness_analysis_pathways.jl` and
`static_options.jl`. Each function builds a Figure and *saves* it to
`pd_config["plot_output_path"]` (px_per_unit=2), then `@info`s the filename.

This file is *included by* the main scripts and assumes `include("src/common.jl")` and a Makie
backend (`using CairoMakie`) are already in scope. Functions read presentation styling passed as
arguments plus a few analysis globals defined by the including script (`n_options`, and, for the
static plots, `scenario_names`, `dhw_scenarios`, `n_dhw`, `option_colors`, `option_labels`,
`dhw_linestyles`, `dhw_markers`, `dhw_model_names`).
=#

# ── Shared presentation constants ─────────────────────────────────────────────

# Pathways figures (9 = GBR/worst/best × 3 metrics, grouped by metric; 3 = metric panels)
metric_labels = [
    "Years >20% coral cover\n(GBR)", "Years >20% coral cover\n(worst reefs)", "Years >20% coral cover\n(best reefs)",
    "Cumulative cover\n(GBR)", "Cumulative cover\n(worst reefs)", "Cumulative cover\n(best reefs)",
    "Cumulative evenness\n(GBR)", "Cumulative evenness\n(worst reefs)", "Cumulative evenness\n(best reefs)"
]
metric_titles = ["Years >20% coral cover", "Cumulative cover", "Cumulative evenness"]

# Static-option figures
aggreg_metric_labels = [
    "Years >20%\n(worst)", "Years >20%\n(best)", "Years >20%\n(net)",
    "Cum. cover\n(worst)", "Cum. cover\n(best)", "Cum. cover\n(net)",
    "Cum. evenness\n(worst)", "Cum. evenness\n(best)", "Cum. evenness\n(net)"
]
boxplot_metric_labels = ["Years >20% coral cover", "Cumulative cover", "Cumulative evenness"]

# ── Parameter-set filename / label helpers ────────────────────────────────────

# Compact scientific notation for seed budgets, e.g. 1_000_000 → "1e6".
_sci(n) = (e = floor(Int, log10(n)); c = n / 10^e; isinteger(c) ? "$(round(Int,c))e$e" : "$(round(c; digits=1))e$e")
# Filename suffix and human title for a parameter set NamedTuple `(dhw, N_seed, n_locations)`.
_param_suffix(ps) = "dhw$(ps.dhw)_n$(_sci(ps.N_seed))_l$(ps.n_locations)"
_param_title(ps) = "dhw $(ps.dhw), N_seed $(_sci(ps.N_seed)), $(ps.n_locations) locs"

# ── Shared drawing helpers ────────────────────────────────────────────────────

# Horizontal dodge offset for option `o_i` of `n` options sharing a group of total width `width`.
_dodge_offset(o_i, n, width) = (o_i - (n + 1) / 2) * (width / n)

"""
    _draw_option_series!(ax, x, point, p10, p90; color, marker=:circle, markersize=10,
        whiskerwidth=6, line=false)

Draw one starting option's point-with-whisker series: markers at `(x, point)` with P10–P90
errorbars. `line=true` links the markers (`scatterlines!`) for across-parameter trends; callers
using it must pass `x` already sorted. Shared by Figures B, C, the B/C comparison and the
parameter scatter so the dodge/marker/errorbar styling stays in one place.
"""
function _draw_option_series!(
    ax, x, point, p10, p90;
    color, marker=:circle, markersize=10, whiskerwidth=6, line=false
)
    if line
        scatterlines!(ax, x, point; color, marker, markersize)
    else
        scatter!(ax, x, point; color, marker, markersize)
    end
    errorbars!(ax, x, point, point .- p10, p90 .- point; color, whiskerwidth)
end

# Single-section "Starting option" colour legend, shared by Figures B, C and the parameter scatter.
function _option_legend!(pos, option_colors, option_labels)
    Legend(pos,
        [MarkerElement(; marker=:circle, color=option_colors[i]) for i in eachindex(option_labels)],
        option_labels;
        title="Starting option", framevisible=false
    )
end

# ── Pathways figures ──────────────────────────────────────────────────────────

"""
    plot_pathways_boxplot(option_names, option_pathways, opt_metric_mats, cf_metric_full,
        option_colors, option_labels, ps)

Figure A — per starting option, boxplots of per-reef deltas (opt − cf) pooled over all reefs ×
all downstream pathways, one panel per base metric. `ps` is the parameter-set NamedTuple
`(dhw, N_seed, n_locations)`.
"""
function plot_pathways_boxplot(
    option_names, option_pathways, opt_metric_mats, cf_metric_full,
    option_colors, option_labels, ps
)
    fig = Figure(size=(1200, 450))
    for (m_i, title) in enumerate(metric_titles)
        ax = Axis(fig[1, m_i];
            xticks=(1:n_options, option_labels),
            xticklabelrotation=π / 4,
            ylabel=m_i == 1 ? "Difference in performance vs counterfactual (per reef)" : "",
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
    fname = "robustness_pathways_boxplot_$(_param_suffix(ps)).png"
    save(joinpath(pd_config["plot_output_path"], fname), fig; px_per_unit=2)
    @info "Saved $(fname)"
end

"""
    plot_pathways_cvar(cvar_df, option_names, option_colors, option_labels, ps)

Figure B — per starting option, median ratio with P10–P90 whiskers across the 9 metric variants
(GBR/worst/best per metric, grouped by metric), dodged by option. `cvar_df` is the tidy table
from [`pathways_cvar_ranges`]; `ps` is the parameter-set NamedTuple `(dhw, N_seed, n_locations)`.
"""
function plot_pathways_cvar(
    cvar_df, option_names, option_colors, option_labels, ps
)
    fig = Figure(size=(1300, 480))
    ax = Axis(fig[1, 1];
        xticks=(1:9, metric_labels),
        xticklabelsize=11,
        xticklabelrotation=π / 4,
        ylabel="Relative performance against counterfactual"
    )
    hlines!(ax, [0]; color=:black, linewidth=2)
    vlines!(ax, [3.5, 6.5]; color=:gray70, linewidth=1, linestyle=:dash)

    for (o_i, option) in enumerate(option_names)
        df = cvar_df[cvar_df.start_option .== string(option), :]
        isempty(df) && continue
        sort!(df, :metric_idx)

        x = df.metric_idx .+ _dodge_offset(o_i, n_options, 0.6)
        _draw_option_series!(ax, x, df.median, df.p10, df.p90; color=option_colors[o_i])
    end

    _option_legend!(fig[1, 2], option_colors, option_labels)

    fname = "robustness_pathways_$(_param_suffix(ps)).png"
    save(joinpath(pd_config["plot_output_path"], fname), fig; px_per_unit=2)
    @info "Saved $(fname)"
end

"""
    plot_pathways_weighted(weighted_tail_stats, option_names, option_colors, option_labels,
        ps; point_stat=:median)

Figure C — same layout as Figure B, but points come from the probability-weighted distribution:
point = `point_stat` (:median or :mean), whiskers span P10–P90. `ps` is the parameter-set
NamedTuple `(dhw, N_seed, n_locations)`.
"""
function plot_pathways_weighted(
    weighted_tail_stats, option_names, option_colors, option_labels, ps;
    point_stat::Symbol=:median
)
    fig = Figure(size=(1300, 480))
    ax = Axis(fig[1, 1];
        xticks=(1:9, metric_labels),
        xticklabelsize=11,
        xticklabelrotation=π / 4,
        ylabel="Probability weighted performance\nagainst counterfactual",
    )
    hlines!(ax, [0]; color=:black, linewidth=2)
    vlines!(ax, [3.5, 6.5]; color=:gray70, linewidth=1, linestyle=:dash)

    for (o_i, option) in enumerate(option_names)
        df = weighted_tail_stats[weighted_tail_stats.start_option .== string(option), :]
        isempty(df) && continue
        sort!(df, :metric_idx)

        x = df.metric_idx .+ _dodge_offset(o_i, n_options, 0.6)
        _draw_option_series!(ax, x, df[!, point_stat], df.p10, df.p90; color=option_colors[o_i])
    end

    _option_legend!(fig[1, 2], option_colors, option_labels)

    fname = "robustness_pathways_weighted_$(point_stat)_$(_param_suffix(ps)).png"
    save(joinpath(pd_config["plot_output_path"], fname), fig; px_per_unit=2)
    @info "Saved $(fname)"
end

"""
    plot_robustness_param_scatter(rob_df, option_names, option_colors, option_labels)

Worst-over-DHW pathway robustness per (starting option × parameter set). x = parameter set,
y = robustness, colour = starting option; point = median over pathways, whiskers = P10/P90, one
line linking each starting option across parameter sets. `rob_df` is the output of
[`worst_dhw_robustness`].
"""
function plot_robustness_param_scatter(rob_df, option_names, option_colors, option_labels)
    # Parameter-set x-axis: unique (N_seed, n_locations), sorted by N_seed then n_locations.
    combos = sort(unique([(r.N_seed, r.n_locations) for r in eachrow(rob_df)]))
    combo_idx = Dict(c => i for (i, c) in enumerate(combos))
    combo_labels = ["$(_sci(c[1])) seeds · $(c[2]) locs" for c in combos]

    fig = Figure(size=(900, 480))
    ax = Axis(fig[1, 1];
        xticks=(1:length(combos), combo_labels),
        xlabel="Parameter set",
        ylabel="Worst-case performance",
        xticklabelrotation = π/4
    )
    hlines!(ax, [0]; color=:black, linewidth=2)

    for (o_i, option) in enumerate(option_names)
        df = rob_df[rob_df.start_option .== string(option), :]
        isempty(df) && continue
        x = [combo_idx[(r.N_seed, r.n_locations)] for r in eachrow(df)] .+
            _dodge_offset(o_i, n_options, 0.5)
        order = sortperm(x)
        _draw_option_series!(ax, x[order], df.median[order], df.p10[order], df.p90[order];
            color=option_colors[o_i], line=true)
    end

    _option_legend!(fig[1, 2], option_colors, option_labels)

    save(
        joinpath(
            pd_config["plot_output_path"],
            "robustness_pathways_param_scatter.png"
        ),
        fig; px_per_unit=2
    )
    @info "Saved robustness_pathways_param_scatter.png"
end

"""
    plot_robustness_vs_diversity(rob_div_df, option_names, option_colors, option_labels)

Facetted scatter of robustness (y) vs pathway diversity (x), coloured by starting option, one
panel per parameter set. Point = median robustness over pathways with P10–P90 whiskers (same
pattern as [`plot_robustness_param_scatter`]); dominated (non-Pareto-optimal) options are hollow.
Panels are gridded so columns = number of locations and rows = number of seeds. All panels share
the same x scale and the same y scale (limits fixed globally, and wide enough to contain the
whiskers). Each panel is annotated with a Kendall τ that folds in the P10/median/P90 robustness
spread via [`panel_tau_diversity_robustness`]. `rob_div_df` is the output of
[`join_robustness_diversity`].
"""
function plot_robustness_vs_diversity(
    rob_div_df, option_names, option_colors, option_labels
)
    n_seeds = sort(unique(rob_div_df.N_seed))       # rows (grid y-axis = number of seeds)
    n_locs = sort(unique(rob_div_df.n_locations))   # columns (grid x-axis = number of locations)
    n_r, n_c = length(n_seeds), length(n_locs)

    # Shared limits: same x scale and same y scale across every panel (with a small margin).
    _pad(lo, hi) = (m = 0.05 * (hi - lo + eps()); (lo - m, hi + m))
    xlims = _pad(extrema(rob_div_df.pathway_diversity)...)
    # y range must span the whiskers, not just the median points, so they are never clipped.
    ylims = _pad(minimum(rob_div_df.robustness_p10), maximum(rob_div_df.robustness_p90))

    fig = Figure(size=(320 * n_c + 220, 260 * n_r + 90))
    for (ri, ns) in enumerate(n_seeds)
        Label(fig[ri, 0], "$(_sci(ns)) seeds"; rotation=π / 2, font=:bold, tellheight=false)
    end
    for (ci, nl) in enumerate(n_locs)
        Label(fig[0, ci], "$(nl) locations"; font=:bold, tellwidth=false)
    end

    for (ri, ns) in enumerate(n_seeds), (ci, nl) in enumerate(n_locs)
        if ns == 1e6
            ylims = (-0.0004, 0.002)
        else
            ylims = _pad(minimum(rob_div_df.robustness_p10), maximum(rob_div_df.robustness_p90))
        end
        ax = Axis(fig[ri, ci];
            limits=(xlims, ylims),
            xlabel=ri == n_r ? "Worst-case pathway diversity" : "",
            ylabel=ci == 1 ? "Worst-case performance" : ""
        )
        sub = rob_div_df[(rob_div_df.N_seed .== ns) .& (rob_div_df.n_locations .== nl), :]
        for (o_i, option) in enumerate(option_names)
            r = sub[sub.start_option .== string(option), :]
            isempty(r) && continue
            # Dominated (non-Pareto-optimal) options are drawn as hollow circles.
            dom = r[!, "dominated?"]
            nd = .!dom
            any(nd) && scatter!(ax, r.pathway_diversity[nd], r.robustness[nd];
                color=option_colors[o_i], markersize=12)
            any(dom) && scatter!(ax, r.pathway_diversity[dom], r.robustness[dom];
                color=:transparent, strokecolor=option_colors[o_i], strokewidth=1.5, markersize=12)
            errorbars!(
                ax, r.pathway_diversity, r.robustness,
                r.robustness .- r.robustness_p10, r.robustness_p90 .- r.robustness;
                color=option_colors[o_i], whiskerwidth=6
            )
        end
        # Descriptive rank agreement for this panel, folding in the P10/median/P90 robustness spread.
        tau = panel_tau_diversity_robustness(sub)
        text!(ax, 0.03, 0.97;
            text=isnan(tau) ? "Kendall-τ = n/a" : "Kendall-τ = $(round(tau; digits=2))",
            space=:relative, align=(:left, :top), fontsize=12
        )
    end

    Legend(fig[1:n_r, n_c + 1],
        [
            [MarkerElement(; marker=:circle, color=option_colors[i]) for i in 1:n_options],
            [
                MarkerElement(; marker=:circle, color=:gray40),
                MarkerElement(; marker=:circle, color=:transparent, strokecolor=:gray40, strokewidth=1.5)
            ]
        ],
        [option_labels, ["Dominating", "Dominated"]],
        ["Starting option", "Marker"];
        framevisible=false
    )

    save(
        joinpath(
            pd_config["plot_output_path"],
            "robustness_vs_diversity.png"
        ),
        fig; px_per_unit=2
    )
    @info "Saved robustness_vs_diversity.png"
end

# Figure B/C comparison: worst/best variants only (GBR dropped), grouped by metric.
# `comparison_metric_idxs` are positions in the 9-variant `cvar_metric_labels` vocabulary.
comparison_metric_idxs = [2, 3, 5, 6, 8, 9]
comparison_metric_labels = [
    "Years >20%\n(worst)", "Years >20%\n(best)",
    "Cum. cover\n(worst)", "Cum. cover\n(best)",
    "Cum. evenness\n(worst)", "Cum. evenness\n(best)"
]

"""
    plot_pathways_weighting_comparison(cvar_df, weighted_tail_stats, option_names,
        option_colors, option_labels, ps; point_stat=:median)

Figures B and C overlaid — for each starting option × metric variant, the **unweighted**
per-pathway range (Figure B, circle) and the **probability-weighted** range (Figure C, diamond)
are drawn side by side with a small horizontal offset, connected by a thin line, so their
absolute values *and* the shift from probability-weighting are both readable. Both series keep
their P10–P90 whiskers. The GBR variant is dropped: only worst/best × 3 metrics (6 x positions).
`point_stat` (:median or :mean) selects the weighted point statistic; the unweighted point is
always the median. `ps` is the parameter-set NamedTuple `(dhw, N_seed, n_locations)`.
"""
function plot_pathways_weighting_comparison(
    cvar_df, weighted_tail_stats, option_names, option_colors, option_labels, ps;
    point_stat::Symbol=:median
)
    # Map the 6 kept metric variants to compact x positions 1..6.
    xpos = Dict(mi => p for (p, mi) in enumerate(comparison_metric_idxs))

    fig = Figure(size=(1300, 520))
    ax = Axis(fig[1, 1];
        xticks=(1:length(comparison_metric_idxs), comparison_metric_labels),
        xticklabelsize=11,
        xticklabelrotation=π / 4,
        ylabel="Relative performance against counterfactual"
    )
    hlines!(ax, [0]; color=:black, linewidth=2)
    vlines!(ax, [2.5, 4.5]; color=:gray70, linewidth=1, linestyle=:dash)

    pair_offset = 0.055  # half-gap between the unweighted and weighted markers of one option

    for (o_i, option) in enumerate(option_names)
        color = option_colors[o_i]
        offset = _dodge_offset(o_i, n_options, 0.72)

        uw = cvar_df[cvar_df.start_option .== string(option), :]
        wt = weighted_tail_stats[weighted_tail_stats.start_option .== string(option), :]
        (isempty(uw) || isempty(wt)) && continue
        uw = uw[in.(uw.metric_idx, Ref(comparison_metric_idxs)), :]
        wt = wt[in.(wt.metric_idx, Ref(comparison_metric_idxs)), :]
        sort!(uw, :metric_idx)
        sort!(wt, :metric_idx)

        base = [xpos[mi] for mi in uw.metric_idx] .+ offset
        x_uw = base .- pair_offset
        x_wt = [xpos[mi] for mi in wt.metric_idx] .+ offset .+ pair_offset

        # Connector between the two estimates of the same option × metric variant.
        for k in eachindex(base)
            lines!(ax, [x_uw[k], x_wt[k]], [uw.median[k], wt[!, point_stat][k]];
                color=(color, 0.5), linewidth=1)
        end

        # Unweighted (Figure B) — filled circle. Probability-weighted (Figure C) — filled diamond.
        _draw_option_series!(ax, x_uw, uw.median, uw.p10, uw.p90;
            color, marker=:circle, markersize=10, whiskerwidth=5)
        _draw_option_series!(ax, x_wt, wt[!, point_stat], wt.p10, wt.p90;
            color, marker=:diamond, markersize=11, whiskerwidth=5)
    end

    Legend(fig[1, 2],
        [
            [MarkerElement(; marker=:circle, color=option_colors[i]) for i in 1:n_options],
            [MarkerElement(; marker=:circle, color=:gray40),
             MarkerElement(; marker=:diamond, color=:gray40)]
        ],
        [option_labels, ["Unweighted", "Prob. weighted ($(point_stat))"]],
        ["Starting option", "Weighting"];
        framevisible=false
    )

    fname = "robustness_pathways_comparison_$(point_stat)_$(_param_suffix(ps)).png"
    save(joinpath(pd_config["plot_output_path"], fname), fig; px_per_unit=2)
    @info "Saved $(fname)"
end

# ── Static-option figures ─────────────────────────────────────────────────────

"""
    plot_robustness(aggreg, N_seed_total, min_iv, rcp, horizon)

Aggregate 9-column line figure from `cvar_aggregation`'s `aggreg` YAXArray.
"""
function plot_robustness(aggreg, N_seed_total, min_iv, rcp, horizon)
    nseed_label = "1e$(round(Int, log10(N_seed_total)))"

    fig = Figure(size=(1300, 480))
    ax = Axis(fig[1, 1];
        xticks=(1:9, aggreg_metric_labels),
        xticklabelsize=10,
        ylabel="Relative performance against \ncounterfactual (no interv.)",
        title="Relative performance Δ/counterfactual - RCP $(rcp), N_seed $(nseed_label), min_locs $(min_iv), horizon $(horizon) yrs"
    )

    hlines!(ax, [0]; color=:black, linewidth=2)
    vlines!(ax, [3.5, 6.5]; color=:gray70, linewidth=1, linestyle=:dash)

    for (o_i, _) in enumerate(scenario_names)
        color = option_colors[o_i]
        for (d_i, _) in enumerate(dhw_scenarios)
            y = aggreg.data[o_i, d_i, :]
            lines!(ax, 1:9, y; color, linestyle=dhw_linestyles[d_i], linewidth=1.5)
            scatter!(ax, 1:9, y; color, marker=dhw_markers[d_i], markersize=10)
        end
    end

    Legend(fig[1, 2],
        [LineElement(; color=option_colors[i], linewidth=2) for i in 1:n_options],
        string.(scenario_names);
        title="Option", framevisible=false
    )
    Legend(fig[2, 2],
        [
            [LineElement(; linestyle=dhw_linestyles[i], color=:gray40, linewidth=2),
             MarkerElement(; marker=dhw_markers[i], color=:gray40)]
            for i in 1:n_dhw
        ],
        [dhw_model_names[i] for i in 1:n_dhw];
        title="Global Climate Model", framevisible=false
    )
    rowsize!(fig.layout, 2, Auto(0.4))

    save(
        joinpath(pd_config["plot_output_path"], "robustness_rcp$(rcp)_n$(nseed_label)_m$(min_iv)_h$(horizon).png"),
        fig; px_per_unit=2
    )
    @info "Saved robustness_rcp$(rcp)_n$(nseed_label)_m$(min_iv)_h$(horizon).png"
end

"""
    plot_boxplot(per_reef, N_seed_total, min_iv, rcp, horizon)

Per-reef boxplot (3 panels) from `cvar_aggregation`'s `per_reef` YAXArray, pooling all
reefs across DHW scenarios.
"""
function plot_boxplot(per_reef, N_seed_total, min_iv, rcp, horizon)
    nseed_label = "1e$(round(Int, log10(N_seed_total)))"

    fig = Figure(size=(1200, 450))
    for (m_i, metric) in enumerate(boxplot_metric_labels)
        ax = Axis(fig[1, m_i];
            xticks=(1:n_options, option_labels),
            xticklabelrotation=π / 4,
            ylabel=m_i == 1 ? "Relative difference vs counterfactual" : "",
            title="$(metric)"
        )
        hlines!(ax, [0]; color=:black, linewidth=2)

        for (o_i, _) in enumerate(scenario_names)
            # pool all reefs × all dhw_scenarios for this option
            y = vec(per_reef[option=o_i, metric=m_i])
            boxplot!(ax, fill(o_i, length(y)), y; color=option_colors[o_i])
        end
    end

    save(
        joinpath(pd_config["plot_output_path"], "robustness_boxplot_rcp$(rcp)_n$(nseed_label)_m$(min_iv)_h$(horizon).png"),
        fig; px_per_unit=2
    )
    @info "Saved robustness_boxplot_rcp$(rcp)_n$(nseed_label)_m$(min_iv)_h$(horizon).png"
end

"""
    plot_net_summary(summary_stats, config_labels, net_syms, net_labels, net_colors, horizon)

Net-metric summary: P10 / median / P90 of each net metric over all options × dhw_scenario, one
point per parameter set (config), whiskers = P10–P90.
"""
function plot_net_summary(summary_stats, config_labels, net_syms, net_labels, net_colors, horizon)
    n_cfg = length(config_labels)

    fig = Figure(size=(1100, 500))
    ax = Axis(fig[1, 1];
        xticks=(1:n_cfg, config_labels),
        xticklabelsize=10,
        ylabel="Net relative performance vs counterfactual",
        title="Net metrics (P10 / median / P90 over options × DHW) — horizon $(horizon) yrs"
    )
    hlines!(ax, [0]; color=:black, linewidth=1)

    dodge = length(net_syms) == 1 ? [0.0] : collect(range(-0.2, 0.2; length=length(net_syms)))
    for (m_i, m) in enumerate(net_syms)
        x = (1:n_cfg) .+ dodge[m_i]
        med = summary_stats[metric=ADRIA.At(m), stat=ADRIA.At(:median)].data
        lo = med .- summary_stats[metric=ADRIA.At(m), stat=ADRIA.At(:p10)].data
        hi = summary_stats[metric=ADRIA.At(m), stat=ADRIA.At(:p90)].data .- med
        errorbars!(ax, x, med, lo, hi; color=net_colors[m_i], whiskerwidth=8)
        scatter!(ax, x, med; color=net_colors[m_i], markersize=10)
    end

    Legend(fig[1, 2],
        [MarkerElement(; marker=:circle, color=net_colors[i]) for i in 1:length(net_syms)],
        net_labels; title="Net metric", framevisible=false
    )

    save(
        joinpath(pd_config["plot_output_path"], "robustness_net_summary_h$(horizon).png"),
        fig; px_per_unit=2
    )
    @info "Saved robustness_net_summary_h$(horizon).png"
end
