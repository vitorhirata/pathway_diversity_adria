#=
Visualisation for the static-option detail plots: GBR choropleth maps (total seeds,
Δ years-above-20%-cover, Δ cumulative cover), seeding-frequency histograms, per-option
time-series, and the per-scenario seeding-map GIF.

Each plotting function builds a Figure and *saves* it to `pd_config["plot_output_path"]`
(px_per_unit=2) with a parameter-set filename suffix, then `@info`s the filename. Functions take
the per-parameter-set data (`sel` from `select_static_scenarios`, and the derived arrays) as
arguments, and read stable styling/label globals defined by the including script
(`scenario_names`, `dhw_model_names`, `all_centroids`).

This file is *included by* the main script and assumes `include("src/common.jl")`, a Makie backend
(`using CairoMakie`), `GeoMakie`, and `NaturalEarth` are already in scope. It also reuses `_sci`
from `visualization/robustness.jl`, which must be included first. Map style follows
`GeoMakie Map Style Guide.md`.
=#

# ── Parameter-set filename suffix ─────────────────────────────────────────────
# e.g. (RCP="45", N_seed=1e6, min_iv=200) → "rcp45_n1e6_l200"
_static_suffix(cfg) = "rcp$(cfg.RCP)_n$(_sci(cfg.N_seed))_l$(cfg.min_iv)"

# ── Shared GBR map infrastructure ─────────────────────────────────────────────

# Pre-load NaturalEarth datasets (cached to disk after first download)
_ne_land = naturalearth("land", 10)
_ne_places = naturalearth("populated_places", 10)

# Shared GBR map extent
_gbr_lon_min, _gbr_lon_max = 141.8, 153.7
_gbr_lat_min, _gbr_lat_max = -25.2, -9.8

function _gbr_annotations!(ax)
    # Step 3: city labels
    for feat in _ne_places
        p = feat.properties
        lon = get(p, :LONGITUDE, nothing)
        lat = get(p, :LATITUDE, nothing)
        (isnothing(lon) || isnothing(lat)) && continue
        _gbr_lon_min <= lon <= _gbr_lon_max || continue
        _gbr_lat_min <= lat <= _gbr_lat_max || continue
        get(p, :ADM0NAME, "") == "Australia" || continue
        get(p, :SCALERANK, 99) <= 6 || continue
        scatter!(ax, [lon], [lat]; color=:black, markersize=5)
        text!(ax, lon - 0.08, lat;
            text=get(p, :NAME, ""), fontsize=8, align=(:right, :center), color=:gray20)
    end

    # Step 4: scale bar
    bar_lat = _gbr_lat_min + 0.5
    bar_lon0 = _gbr_lon_min + 0.3
    bar_lon1 = bar_lon0 + 100.0 / (111.32 * cosd(abs(bar_lat)))
    cap_h = 0.07
    lines!(ax, [bar_lon0, bar_lon1], [bar_lat, bar_lat]; color=:black, linewidth=2.5)
    lines!(
        ax,
        [bar_lon0, bar_lon0],
        [bar_lat - cap_h, bar_lat + cap_h];
        color=:black,
        linewidth=2.5
    )
    lines!(
        ax,
        [bar_lon1, bar_lon1],
        [bar_lat - cap_h, bar_lat + cap_h];
        color=:black,
        linewidth=2.5
    )
    text!(
        ax,
        bar_lon0,
        bar_lat - cap_h - 0.08;
        text="0",
        align=(:center, :top),
        fontsize=9,
        color=:black
    )
    text!(
        ax,
        bar_lon1,
        bar_lat - cap_h - 0.08;
        text="100 km",
        align=(:center, :top),
        fontsize=9,
        color=:black
    )

    # Step 5: north arrow
    arr_lon = _gbr_lon_max - 0.6
    arr_lat0 = _gbr_lat_max - 1.5
    arr_dlat = 0.8
    # Makie ≥ 0.24 split `arrows!` into `arrows2d!`/`arrows3d!` and replaced the tip kwargs.
    # SankeyMakie (Makie = "0.21, 0.22") holds this env on 0.22, so support both.
    if isdefined(Makie, :arrows2d!)
        Makie.arrows2d!(ax, [arr_lon], [arr_lat0], [0.0], [arr_dlat];
            color=:black, tiplength=10.32, tipwidth=8.9)
    else
        Makie.arrows!(ax, [arr_lon], [arr_lat0], [0.0], [arr_dlat];
            color=:black, arrowsize=10)
    end
    return text!(ax, arr_lon, arr_lat0 + arr_dlat + 0.12;
        text="N", align=(:center, :bottom), fontsize=13, font=:bold)
end

"""
    panel_map_figure(panel_data, n_panels; title_fn, cbar_label, colormap, colorrange,
                     gray_points=nothing)

GBR scatter map of `n_panels` panels side by side sharing a single colorbar.
`panel_data(p)` returns `(points, values)` for panel `p`; the optional `gray_points(p)`
returns locations drawn as inactive background.
"""
function panel_map_figure(
    panel_data,
    n_panels;
    title_fn,
    cbar_label,
    colormap,
    colorrange,
    gray_points=nothing
)
    fig = Figure(; size=(600 * n_panels, 600))
    for p in 1:n_panels
        ax = GeoAxis(
            fig[1, p];
            dest="+proj=longlat +datum=WGS84",
            limits=(_gbr_lon_min, _gbr_lon_max, _gbr_lat_min, _gbr_lat_max),
            title=title_fn(p),
            xgridcolor=(:gray, 0.15),
            ygridcolor=(:gray, 0.15),
            xticklabelsize=7,
            yticklabelsize=7
        )
        poly!(
            ax,
            _ne_land.geometry;
            color=RGBf(0.93, 0.91, 0.87),
            strokewidth=0.5,
            strokecolor=:gray40
        )
        if !isnothing(gray_points)
            scatter!(ax, gray_points(p); color=:gray80, markersize=4, alpha=0.5)
        end

        points, vals = panel_data(p)
        order = sortperm(vals)
        scatter!(ax, points[order]; color=vals[order], colormap=colormap,
            colorrange=colorrange, markersize=4, alpha=0.7)

        _gbr_annotations!(ax)
    end
    Colorbar(fig[1, n_panels + 1];
        colorrange=colorrange,
        colormap=colormap,
        label=cbar_label,
        height=Relative(0.65)
    )
    rowsize!(fig.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig)

    return fig
end

# ── Seeding-frequency histograms ──────────────────────────────────────────────

"""Seeding-frequency histograms per intervention scenario (one row per DHW member)."""
function plot_seeding_frequency_per_option(seed_per_reef_per_ts_scen, sel, cfg)
    fig = Figure(; size=(length(scenario_names) * 350, sel.n_dhw * 400))
    for d in 1:sel.n_dhw, (opt_i, scen_name) in enumerate(scenario_names)
        seeded_binary = seed_per_reef_per_ts_scen[scenarios=sel.iv_col(d, opt_i)] .> 0
        seeding_freq = vec(mean(seeded_binary; dims=1)) .* 100
        ax = Axis(
            fig[d, opt_i];
            xlabel="Seeding frequency (%)",
            ylabel=opt_i == 1 ? "Number of reefs\n$(dhw_model_names[d])" : "",
            title="$scen_name — $(dhw_model_names[d])",
            xticks=0:20:100
        )
        hist!(ax, seeding_freq; bins=0:5:100)
    end
    fname = "seeding_frequency_per_option_$(_static_suffix(cfg)).png"
    save(joinpath(pd_config["plot_output_path"], fname), fig)
    @info "Saved $(fname)"
end

"""Seeding-frequency histogram aggregating all intervention scenarios and DHW members."""
function plot_seeding_frequency_all(seed_per_reef_per_ts_scen, cfg)
    seeding_freq_all = vec(mean(seed_per_reef_per_ts_scen.data .> 0; dims=(1, 3))) .* 100
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        xlabel="Seeding frequency (% of timestep–scenario combinations)",
        ylabel="Number of reefs",
        title="Seeding frequency — all intervention scenarios, all DHW members",
        xticks=0:10:100
    )
    hist!(ax, seeding_freq_all; bins=0:5:100)
    fname = "seeding_frequency_all_$(_static_suffix(cfg)).png"
    save(joinpath(pd_config["plot_output_path"], fname), fig)
    @info "Saved $(fname)"
end

# ── Total seeds per location ──────────────────────────────────────────────────

"""
    plot_total_seeds(total_seeds, sel, all_centroids; seeds_dhw, n_seed_cols, colorrange,
                     colgap, rowgap, panel_height)

Grid of GBR maps of total deployed coral per location, one panel per option (title = option
name), all sharing a single colorbar, for a single DHW member `seeds_dhw`.

Panel cells are sized to the GBR map's own aspect ratio so there is no whitespace baked into
each axis; `colgap`/`rowgap` (in px) then set the actual spacing between panels.
"""
function plot_total_seeds(
    total_seeds,
    sel,
    all_centroids;
    seeds_dhw=1,
    n_seed_cols=3,
    colorrange=(0.0, Float64(maximum(total_seeds))),
    colgap=6,
    rowgap=6,
    panel_height=380
)
    n_rows = cld(length(scenario_names), n_seed_cols)
    fig = Figure()
    for (opt_i, scen_name) in enumerate(scenario_names)
        row = cld(opt_i, n_seed_cols)
        col = mod1(opt_i, n_seed_cols)
        ax = GeoAxis(
            fig[row, col];
            dest="+proj=longlat +datum=WGS84",
            limits=(_gbr_lon_min, _gbr_lon_max, _gbr_lat_min, _gbr_lat_max),
            title=uppercasefirst(replace(string(scen_name), '_' => ' ')),
            xgridcolor=(:gray, 0.15),
            ygridcolor=(:gray, 0.15),
            xticklabelsize=7,
            yticklabelsize=7
        )
        # Keep y ticks only on the left column and x ticks only on the bottom row
        col == 1 || hideydecorations!(ax; grid=false)
        row == n_rows || hidexdecorations!(ax; grid=false)
        poly!(
            ax,
            _ne_land.geometry;
            color=RGBf(0.93, 0.91, 0.87),
            strokewidth=0.5,
            strokecolor=:gray40
        )
        scatter!(
            ax,
            all_centroids[total_seeds[:, sel.iv_col(seeds_dhw, opt_i)] .== 0];
            color=:gray80, markersize=4, alpha=0.5
        )

        has_seeds = total_seeds[:, sel.iv_col(seeds_dhw, opt_i)] .> 0
        vals = total_seeds[has_seeds, sel.iv_col(seeds_dhw, opt_i)]
        order = sortperm(vals)
        scatter!(ax, all_centroids[has_seeds][order]; color=vals[order], colormap=:viridis,
            colorrange=colorrange, markersize=4, alpha=0.7)

        _gbr_annotations!(ax)
    end
    Colorbar(fig[:, n_seed_cols + 1];
        colorrange=colorrange,
        colormap=:viridis,
        label="Total deployed coral",
        height=Relative(0.65)
    )

    # Match each map panel to the GBR extent's aspect (cos-corrected for latitude) so the
    # axis has no internal horizontal whitespace; only then do the gaps control spacing.
    map_aspect = (_gbr_lon_max - _gbr_lon_min) * cosd((_gbr_lat_min + _gbr_lat_max) / 2) /
                 (_gbr_lat_max - _gbr_lat_min)
    for r in 1:n_rows
        rowsize!(fig.layout, r, Fixed(panel_height))
    end
    for c in 1:n_seed_cols
        colsize!(fig.layout, c, Aspect(1, map_aspect))
    end
    colgap!(fig.layout, colgap)
    rowgap!(fig.layout, rowgap)
    resize_to_layout!(fig)

    return fig
end

"""Total seeds per location map (one panel per option) for a single DHW member, saved."""
function plot_total_seeds_map(total_seeds, sel, cfg; seeds_dhw=1, colorrange=(0.0, 1e6))
    fig = plot_total_seeds(total_seeds, sel, all_centroids; seeds_dhw, colorrange)
    fname = "total_seeds_$(_static_suffix(cfg)).png"
    save(joinpath(pd_config["plot_output_path"], fname), fig; px_per_unit=2)
    @info "Saved $(fname)"
end

# ── Performance-difference maps ───────────────────────────────────────────────

# Each performance map is drawn once per location filter: filter 1 keeps only the reefs this
# scenario/DHW actually seeded, filter 2 keeps every reef seeded by at least one scenario (so all
# options share a common reef set, and reefs this option skipped are still shown).
function _plot_perf_maps(
    diff_mat, total_seeds, active_mask, sel, cfg;
    fname_stem, title_prefix, cbar_label, colormap, colorrange
)
    filter_labels = ["seeded reefs", "all reefs"]
    filter_names = ["seeded_reefs", "all_reefs"]
    for (opt_i, scen_name) in enumerate(scenario_names), f in 1:2
        # The seeded-reef filter depends on the DHW member, so resolve it per panel
        mask(d) = f == 1 ? (total_seeds[:, sel.iv_col(d, opt_i)] .> 0) : active_mask

        fig = panel_map_figure(
            d -> (all_centroids[mask(d)], diff_mat[mask(d), sel.iv_col(d, opt_i)]),
            sel.n_dhw;
            title_fn=d -> "$title_prefix — $scen_name vs counterfactual\n$(dhw_model_names[d]) — $(filter_labels[f])",
            cbar_label=cbar_label,
            colormap=colormap,
            colorrange=colorrange
        )
        fname = "$(fname_stem)_$(scen_name)_$(filter_names[f])_$(_static_suffix(cfg)).png"
        save(joinpath(pd_config["plot_output_path"], fname), fig; px_per_unit=2)
        @info "Saved $(fname)"
    end
end

"""Δ years-above-20%-cover maps (per option × location filter), one panel per DHW member."""
function plot_nyrs_above_maps(n_yrs_diff, total_seeds, active_mask, sel, cfg)
    _plot_perf_maps(
        n_yrs_diff, total_seeds, active_mask, sel, cfg;
        fname_stem="n_yrs_above_target",
        title_prefix="Δ years above 20% coral cover",
        cbar_label="Δ years above 20% coral cover (intervention − counterfactual)",
        colormap=Reverse(:vik25),
        colorrange=(-3, 3)
    )
end

"""Δ cumulative-cover maps (per option × location filter), one panel per DHW member."""
function plot_cum_tac_diff_maps(cum_tac_diff, total_seeds, active_mask, sel, cfg)
    # Symmetric colour range from the bulk of the deltas, so a few extreme reefs don't flatten
    # the map (over the active reefs, all intervention columns).
    _cc_lim = quantile(abs.(vec(cum_tac_diff[active_mask, :])), 0.98)
    colorrange = iszero(_cc_lim) ? (-1.0, 1.0) :
        (-round(_cc_lim; sigdigits=2), round(_cc_lim; sigdigits=2))
    _plot_perf_maps(
        cum_tac_diff, total_seeds, active_mask, sel, cfg;
        fname_stem="cum_tac_diff",
        title_prefix="Δ cumulative coral cover",
        cbar_label="Δ cumulative coral cover (km²·years, intervention − counterfactual)",
        colormap=Reverse(:vik),
        colorrange=colorrange
    )
end

# ── Per-option time-series ────────────────────────────────────────────────────

"""
    plot_options_timeseries(rs, metrics, sel, options, cfg)

Per-metric time-series figures (one per metric × DHW member): left panel = each option's scenario
trajectory vs the counterfactual, right panel = option − counterfactual difference.
"""
function plot_options_timeseries(rs, metrics, sel, options, cfg)
    option_names = Symbol.(options.option_name)
    all_names = vcat(option_names, [:unguided_intervention, :no_intervention])
    intervention_names = all_names[1:(end - 1)]
    n_block = length(all_names)         # options + unguided + counterfactual
    n_interv = length(intervention_names)
    # Groups index into a single DHW block (one figure per DHW member)
    scen_groups = Dict{Symbol,BitVector}(
        name => BitVector((1:n_block) .== i) for (i, name) in enumerate(all_names)
    )
    scen_groups_diff = Dict{Symbol,BitVector}(
        name => BitVector((1:n_interv) .== i) for (i, name) in enumerate(intervention_names)
    )

    ts = string.(ADRIA.timesteps(rs))
    tick_pos = collect(1:5:length(ts))
    tick_lbl = ts[1:5:end]
    (length(ts) - 1) % 5 != 0 &&
        (tick_pos = vcat(tick_pos, length(ts)); tick_lbl = vcat(tick_lbl, ts[end]))
    xtick_vals = (tick_pos, tick_lbl)
    xtick_rot = 2 / π

    for (name, metric) in metrics, d in 1:sel.n_dhw
        # This DHW block only, in the order [options..., unguided, counterfactual]
        metric_dhw = ADRIA.DataCube(
            metric.data[:, sel.block_cols(d)];
            timesteps=ADRIA.timesteps(rs),
            scenarios=1:n_block
        )
        metric_diff = ADRIA.DataCube(
            metric_dhw.data[:, 1:(end - 1)] .- metric_dhw.data[:, end];
            timesteps=ADRIA.timesteps(rs),
            scenarios=1:n_interv
        )

        f = Figure(; size=(3200, 800))
        g1 = f[1, 1] = GridLayout()
        g2 = f[1, 2] = GridLayout()

        ax1 = Axis(g1[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot,
            title="$name — $(dhw_model_names[d])")
        ADRIA.viz.scenarios!(g1, ax1, metric_dhw, scen_groups;
            opts=Dict{Symbol,Any}(
                :legend_labels => all_names, :legend => false, :histogram => false
            ))
        ADRIA.viz.scenarios_legend!(g1[1, 0], scen_groups, metric_dhw;
            opts=Dict{Symbol,Any}(:legend_labels => all_names),
            legend_opts=Dict{Symbol,Any}(:padding => (4, 4, 4, 4))
        )

        ax2 = Axis(g2[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot,
            title="$name - counterfactual — $(dhw_model_names[d])",
            ylabel="$name - counterfactual")
        ADRIA.viz.scenarios!(g2, ax2, metric_diff, scen_groups_diff;
            opts=Dict{Symbol,Any}(
                :legend_labels => intervention_names, :legend => false, :histogram => false
            ))

        fname = "options_$(replace(lowercase(name), ' ' => '_'))_dhw$(sel.dhw_scenarios[d])_$(_static_suffix(cfg)).png"
        save(joinpath(pd_config["plot_output_path"], fname), f)
        @info "Saved $(fname)"
    end
end

# ── Selected-locations seeding GIF ────────────────────────────────────────────

"""
    animate_seeding_maps(rs, dom, seed_per_reef_per_ts_scen, sel, cfg; gif_dhw=1)

One animated GBR map per intervention scenario for DHW member `gif_dhw`: seeded reefs (red) over
the seeding window, on the grey reef footprint.
"""
function animate_seeding_maps(rs, dom, seed_per_reef_per_ts_scen, sel, cfg; gif_dhw=1)
    seed_start = Int(rs.inputs.seed_year_start[1])
    n_seed_years = Int(rs.inputs.seed_years[1])
    seed_ts = seed_start:(seed_start + n_seed_years - 1)
    ts_labels = ADRIA.timesteps(rs)[seed_ts]
    plottable_gif = GeoMakie.to_multipoly(dom.loc_data[:, :geometry])

    for (opt_i, scen_name) in enumerate(scenario_names)
        seeded_points = Observable(Point2f[])
        title_obs = Observable("$scen_name — Year: $(ts_labels[1])")

        fig_gif = Figure()
        ga_gif = GeoAxis(
            fig_gif[1, 1];
            dest="+proj=longlat +datum=WGS84",
            limits=(_gbr_lon_min, _gbr_lon_max, _gbr_lat_min, _gbr_lat_max),
            title=title_obs,
            titlesize=20,
            xgridcolor=(:gray, 0.15),
            ygridcolor=(:gray, 0.15),
            xticklabelsize=7,
            yticklabelsize=7,
        )
        poly!(ga_gif, _ne_land.geometry; color=RGBf(0.93, 0.91, 0.87), strokewidth=0.5, strokecolor=:gray40)
        poly!(ga_gif, plottable_gif; color=:gray80)
        scatter!(ga_gif, seeded_points; color=:red, markersize=4)
        _gbr_annotations!(ga_gif)
        rowsize!(fig_gif.layout, 1, Aspect(1, 1.0))
        resize_to_layout!(fig_gif)

        fname = "seeding_map_$(scen_name)_dhw$(sel.dhw_scenarios[gif_dhw])_$(_static_suffix(cfg)).gif"
        record(
            fig_gif,
            joinpath(pd_config["plot_output_path"], fname),
            eachindex(ts_labels);
            framerate=3
        ) do i
            seeded_points[] = all_centroids[seed_per_reef_per_ts_scen[timesteps=i, scenarios=sel.iv_col(gif_dhw, opt_i)] .> 0]
            title_obs[] = "$scen_name — Year: $(ts_labels[i])"
        end
        @info "Saved $(fname)"
    end
end
