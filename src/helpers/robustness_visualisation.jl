#=
Visualisation for the robustness analyses.

Plotting functions for both `robustness_analysis_pathways.jl` and
`robustness_analysis_static_options.jl`. Each function builds a Figure and *saves* it to
`pd_config["plot_output_path"]` (px_per_unit=2), then `@info`s the filename.

This file is *included by* the main scripts and assumes `include("src/common.jl")` and a Makie
backend (`using CairoMakie`) are already in scope. Functions read presentation styling passed as
arguments plus a few analysis globals defined by the including script (`n_options`, and, for the
static plots, `scenario_names`, `dhw_scenarios`, `n_dhw`, `option_colors`, `option_labels`,
`dhw_linestyles`, `dhw_markers`, `dhw_model_names`).
=#

# ── Shared presentation constants ─────────────────────────────────────────────

# Pathways figures (6 = worst/best × 3 metrics; 3 = metric panels)
metric_labels = [
    "Years >20% coral cover\n(worst reefs)", "Years >20% coral cover\n(best reefs)",
    "Cumulative cover\n(worst reefs)", "Cumulative cover\n(best reefs)",
    "Cumulative evenness\n(worst reefs)", "Cumulative evenness\n(best reefs)"
]
metric_titles = ["Years >20% coral cover", "Cumulative cover", "Cumulative evenness"]

# Static-option figures
aggreg_metric_labels = [
    "Years >20%\n(worst)", "Years >20%\n(best)", "Years >20%\n(net)",
    "Cum. cover\n(worst)", "Cum. cover\n(best)", "Cum. cover\n(net)",
    "Cum. evenness\n(worst)", "Cum. evenness\n(best)", "Cum. evenness\n(net)"
]
boxplot_metric_labels = ["Years >20% coral cover", "Cumulative cover", "Cumulative evenness"]

# ── Pathways figures ──────────────────────────────────────────────────────────

"""
    plot_pathways_boxplot(option_names, option_pathways, opt_metric_mats, cf_metric_full,
        option_colors, option_labels, sel_rcp, sel_dhw)

Figure A — per starting option, boxplots of per-reef deltas (opt − cf) pooled over all reefs ×
all downstream pathways, one panel per base metric.
"""
function plot_pathways_boxplot(
    option_names, option_pathways, opt_metric_mats, cf_metric_full,
    option_colors, option_labels, sel_rcp, sel_dhw
)
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
end

"""
    plot_pathways_cvar(med, lo, hi, option_names, option_colors, option_labels, sel_rcp, sel_dhw)

Figure B — per starting option, median CVaR tail ratio with min/max whiskers across the 6 metric
variants (worst/best × 3 metrics), dodged by option.
"""
function plot_pathways_cvar(
    med, lo, hi, option_names, option_colors, option_labels, sel_rcp, sel_dhw
)
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
end

"""
    plot_pathways_weighted(weighted_tail_stats, option_names, option_colors, option_labels,
        sel_rcp, sel_dhw; point_stat=:median)

Figure C — same layout as Figure B, but points come from the probability-weighted distribution:
point = `point_stat` (:median or :mean), whiskers span P10–P90.
"""
function plot_pathways_weighted(
    weighted_tail_stats, option_names, option_colors, option_labels, sel_rcp, sel_dhw;
    point_stat::Symbol=:median
)
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
end

"""
    plot_robustness_param_scatter(rob_df, option_names, option_colors, option_labels, sel_rcp)

Worst-DHW pathway robustness per (starting option × parameter set). x = parameter set,
y = robustness, colour = starting option; point = median over pathways, whiskers = min/max, one
line linking each starting option across parameter sets. `rob_df` is the output of
[`worst_dhw_robustness`].
"""
function plot_robustness_param_scatter(rob_df, option_names, option_colors, option_labels, sel_rcp)
    # Parameter-set x-axis: unique (N_seed, n_locations), sorted by N_seed.
    combos = sort(unique([(r.N_seed, r.n_locations) for r in eachrow(rob_df)]); by=first)
    combo_idx = Dict(c => i for (i, c) in enumerate(combos))
    _sci(n) = (e = floor(Int, log10(n)); c = n / 10^e; isinteger(c) ? "$(round(Int,c))e$e" : "$(round(c; digits=1))e$e")
    combo_labels = ["$(_sci(c[1])) · $(c[2]) locs" for c in combos]

    fig = Figure(size=(900, 480))
    ax = Axis(fig[1, 1];
        xticks=(1:length(combos), combo_labels),
        xlabel="Parameter set",
        ylabel="Worst-DHW robustness\n(median over pathways, min–max)",
        title="Pathway robustness by starting option — RCP $(sel_rcp)"
    )
    hlines!(ax, [0]; color=:black, linewidth=2)

    n_dodge = n_options
    dodge_width = 0.5
    for (o_i, option) in enumerate(option_names)
        df = rob_df[rob_df.start_option .== string(option), :]
        isempty(df) && continue
        color = option_colors[o_i]
        offset = (o_i - (n_dodge + 1) / 2) * (dodge_width / n_dodge)
        x = [combo_idx[(r.N_seed, r.n_locations)] for r in eachrow(df)] .+ offset
        order = sortperm(x)
        x = x[order]
        med = df.median[order]
        scatterlines!(ax, x, med; color, markersize=10)
        errorbars!(
            ax, x, med, med .- df.min[order], df.max[order] .- med;
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
            "robustness_pathways_param_scatter_rcp$(sel_rcp).png"
        ),
        fig; px_per_unit=2
    )
    @info "Saved robustness_pathways_param_scatter_rcp$(sel_rcp).png"
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
