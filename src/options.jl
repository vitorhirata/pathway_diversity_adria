#=
Options experiment: compare five seeding strategy options under two DHW scenario members.

Design:
- 5 options (from option_seed_preference) + unguided + counterfactual = 7 scenarios,
  repeated once per DHW member
- Only dhw_scenario varies; all other parameters are fixed
- Scenario rows are laid out in blocks of `n_scens_per_dhw`, one block per DHW member:
  row (d - 1) * n_scens_per_dhw + s → DHW member d, scenario s.
  The last scenario of each block (s = n_scens_per_dhw) is that block's counterfactual.

To identify the DHW member and option for scenario i after running:
    d = div(i - 1, n_scens_per_dhw) + 1   # dhw_scenarios[d] == rs.inputs.dhw_scenario[i]
    s = mod1(i, n_scens_per_dhw)          # scenario_names[s], or counterfactual if s == 7
=#

include("src/common.jl")
using GeoMakie, GraphMakie, CairoMakie, NaturalEarth

RCP = "45"
seed_years = 20
dhw_scenarios = [7, 10]
n_dhw = length(dhw_scenarios)
n_scens_per_dhw = 7          # 5 options + unguided + counterfactual
n_options = n_scens_per_dhw - 1  # intervention scenarios (options + unguided)

# Row/column helpers. `scen_*` index into the full scenario table (and rs scenario dim);
# `iv_col*` index into intervention-only arrays (counterfactual columns dropped).
scen_rows(d) = ((d - 1) * n_scens_per_dhw + 1):(d * n_scens_per_dhw)
scen_row(d, s) = (d - 1) * n_scens_per_dhw + s
cf_idx(d) = scen_row(d, n_scens_per_dhw)
iv_col(d, o) = (d - 1) * n_options + o

dom = ADRIA.load_domain(pd_config["domain_path"], RCP;
    calib_params_fn=pd_config["coral_param_path"],
    # timeframe: seed_years + 2 (start seeding), 5 (extra years)
    timeframe=(2022, 2022 + seed_years + 2 + 5)
)
fix_common_parameters!(dom)

N_seed_weights = (
    N_seed_TA=0.15, N_seed_CA=0.5, N_seed_CNA=0.0, N_seed_SM=0.35, N_seed_LM=0.0
)
N_seed_total = 1e7

ADRIA.fix_factor!(dom;
    # Environmental
    # Baseline only; every scenario row gets its block's dhw_scenario assigned below.
    dhw_scenario=dhw_scenarios[1],
    # Intervention parameters
    min_iv_locations=200,
    # Seeding parameters
    seed_years=seed_years,
    seeding_devices_per_m2=5,
    a_adapt=5.0,
    N_seed_TA=N_seed_total * N_seed_weights.N_seed_TA,
    N_seed_CA=N_seed_total * N_seed_weights.N_seed_CA,
    N_seed_CNA=N_seed_total * N_seed_weights.N_seed_CNA,
    N_seed_SM=N_seed_total * N_seed_weights.N_seed_SM,
    N_seed_LM=N_seed_total * N_seed_weights.N_seed_LM
)

options = ADRIA.analysis.option_seed_preference(; include_weights=true)
scens_block = []
for option in eachrow(options)
    ADRIA.fix_factor!(dom;
        seed_heat_stress=option[2],
        seed_in_connectivity=option[3],
        seed_out_connectivity=option[4],
        seed_depth=option[5],
        seed_coral_cover=option[6],
        seed_cluster_diversity=option[7],
        seed_geographic_separation=option[8],
        seed_coral_diversity=option[9]
    )
    scen = ADRIA.sample(dom, 2)[1:1, :]
    push!(scens_block, scen)
end

# Unguided seeding
ADRIA.fix_factor!(dom; guided=0)
scen = ADRIA.sample(dom, 2)[1:1, :]
push!(scens_block, scen)

# No seeding
ADRIA.fix_factor!(dom;
    N_seed_TA=0,
    N_seed_CA=0,
    N_seed_CNA=0,
    N_seed_SM=0,
    N_seed_LM=0
)
scen = ADRIA.sample(dom, 2)[1:1, :]
push!(scens_block, scen)

scens_block = vcat(scens_block...)
@assert nrow(scens_block) == n_scens_per_dhw

# One copy of the 7-scenario block per DHW member
scens = repeat(scens_block, n_dhw)
for d in 1:n_dhw
    scens[scen_rows(d), :dhw_scenario] .= dhw_scenarios[d]
end

rs = ADRIA.run_scenarios(dom, scens, RCP)

# Load scenario
# path = "Outputs/"
# rs = ADRIA.load_results(path)

# Reefs that received seeding in every intervention scenario
seed_start = Int(rs.inputs.seed_year_start[1])
n_seed_years = Int(rs.inputs.seed_years[1])
seed_ts = seed_start:(seed_start + n_seed_years - 1)
# Intervention scenarios of every DHW block, ordered so that column iv_col(d, o) holds
# option o under DHW member d.
intervention_scens = vcat([[scen_row(d, o) for o in 1:n_options] for d in 1:n_dhw]...)
seed_per_reef_per_ts_scen = dropdims(
    sum(rs.seed_log[timesteps=seed_ts, scenarios=intervention_scens]; dims=:coral_id);
    dims=:coral_id
)

# always/never seeded: must hold across every timestep AND every scenario (both DHW members)
always_seeded = vec(all(seed_per_reef_per_ts_scen.data .> 0; dims=(1, 3)))
never_seeded = vec(all(seed_per_reef_per_ts_scen.data .== 0; dims=(1, 3)))

# Remove reefs that always or never have seeding
selected_locations = .!(always_seeded .| never_seeded)
selected_locations = findall(selected_locations)

# Validation: seeding coverage across all intervention scenarios
@info("Total reefs  : $(size(dom.loc_data, 1))")
@info("Never seeded : $(sum(never_seeded)) ($(round(100*mean(never_seeded); digits=1))%)")
@info("Always seeded: $(sum(always_seeded)) ($(round(100*mean(always_seeded); digits=1))%)")

s_tac = ADRIA.metrics.scenario_total_cover(rs; locations=selected_locations)
s_rsv = ADRIA.metrics.scenario_rsv(rs; locations=selected_locations)
s_even = ADRIA.metrics.scenario_evenness(rs; locations=selected_locations)

# scenario_relative_juveniles ignores the locations kwarg, so pre-slice manually
_aj = ADRIA.metrics.absolute_juveniles(rs)
_k_area = ADRIA.loc_k_area(rs)[selected_locations]
s_juves = ADRIA.metrics.scenario_relative_juveniles(
    _aj[locations=selected_locations].data, _k_area
)
metrics = Dict(
    "Total absolute cover" => s_tac,
    "Relative Shelter Volume" => s_rsv,
    "Relative Juveniles" => s_juves,
    "Coral Evenness" => s_even
)

# Pre-load NaturalEarth datasets (cached to disk after first download)
_ne_land = naturalearth("land", 10)
_ne_places = naturalearth("populated_places", 10)

# Shared GBR map extent
_gbr_lon_min, _gbr_lon_max = 141.8, 153.7
_gbr_lat_min, _gbr_lat_max = -25.2, -9.8

ts_labels = ADRIA.timesteps(rs)[seed_ts]
all_centroids = ADRIA.centroids(dom.loc_data)
scenario_names = vcat(options.option_name, [:unguided])
dhw_model_names = dom.dhw_scens.properties["model_names"][dhw_scenarios]

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

#=
# Selected locations GIF per intervention scenario, for a single DHW member
gif_dhw = 1
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

    record(
        fig_gif,
        joinpath(pd_config["plot_output_path"], "seeding_map_$(scen_name)_dhw$(dhw_scenarios[gif_dhw]).gif"),
        eachindex(ts_labels);
        framerate=3
    ) do i
        seeded_points[] = all_centroids[seed_per_reef_per_ts_scen[timesteps=i, scenarios=iv_col(gif_dhw, opt_i)] .> 0]
        title_obs[] = "$scen_name — Year: $(ts_labels[i])"
    end
end
=#

# Seeding frequency histograms per intervention scenario (one row per DHW member)
fig_hist = Figure(; size=(length(scenario_names) * 350, n_dhw * 400))
for d in 1:n_dhw, (opt_i, scen_name) in enumerate(scenario_names)
    seeded_binary = seed_per_reef_per_ts_scen[scenarios=iv_col(d, opt_i)] .> 0
    seeding_freq = vec(mean(seeded_binary; dims=1)) .* 100
    ax = Axis(
        fig_hist[d, opt_i];
        xlabel="Seeding frequency (%)",
        ylabel=opt_i == 1 ? "Number of reefs\n$(dhw_model_names[d])" : "",
        title="$scen_name — $(dhw_model_names[d])",
        xticks=0:20:100
    )
    hist!(ax, seeding_freq; bins=0:5:100)
end
save(joinpath(pd_config["plot_output_path"], "seeding_frequency_per_option.png"), fig_hist)

# Seeding frequency histograms aggregating cenarios
seeding_freq_all = vec(mean(seed_per_reef_per_ts_scen.data .> 0; dims=(1, 3))) .* 100
fig_hist_all = Figure()
ax_hist_all = Axis(
    fig_hist_all[1, 1];
    xlabel="Seeding frequency (% of timestep–scenario combinations)",
    ylabel="Number of reefs",
    title="Seeding frequency — all intervention scenarios, all DHW members",
    xticks=0:10:100
)
hist!(ax_hist_all, seeding_freq_all; bins=0:5:100)
save(joinpath(pd_config["plot_output_path"], "seeding_frequency_all.png"), fig_hist_all)

active_mask = .!never_seeded

# Total seeds per location map per intervention scenario (one panel per DHW member)
total_seeds = Array(
    dropdims(
        sum(rs.seed_log[scenarios=intervention_scens]; dims=(:timesteps, :coral_id));
        dims=(:timesteps, :coral_id)
    )
)  # (n_locs, n_dhw * n_options), column iv_col(d, o)
#seed_colorrange = (0.0, Float64(maximum(total_seeds)))
seed_colorrange = (0.0, 1.2e7)

for (opt_i, scen_name) in enumerate(scenario_names)
    fig_seeds = panel_map_figure(
        d -> begin
            has_seeds = total_seeds[:, iv_col(d, opt_i)] .> 0
            (all_centroids[has_seeds], total_seeds[has_seeds, iv_col(d, opt_i)])
        end,
        n_dhw;
        title_fn=d -> "Total deployed coral per location\n$scen_name — $(dhw_model_names[d])",
        cbar_label="Total deployed coral",
        colormap=:viridis,
        colorrange=seed_colorrange,
        gray_points=d -> all_centroids[total_seeds[:, iv_col(d, opt_i)] .== 0]
    )
    save(
        joinpath(pd_config["plot_output_path"], "total_seeds_$(scen_name).png"),
        fig_seeds;
        px_per_unit=2
    )
end

# n_yrs_above_target scatter map per scenario and location filter (one panel per DHW member)
m_tac = Array(ADRIA.metrics.total_absolute_cover(rs)) .* 1e-6  # km²
loc_hab_area_km2 = rs.loc_area .* rs.loc_max_coral_cover .* 1e-6
n_locs_total = size(m_tac, 2)

n_yrs_above = Matrix{Int32}(undef, n_locs_total, nrow(scens))
for s in 1:nrow(scens)
    for l in 1:n_locs_total
        thr = 0.20 * loc_hab_area_km2[l] # Threshold of 20% coral cover
        n_yrs_above[l, s] = Int32(count(m_tac[:, l, s] .>= thr))
    end
end

# Δ vs the counterfactual of the same DHW member
n_yrs_diff = Matrix{Float64}(undef, n_locs_total, nrow(scens))
for d in 1:n_dhw
    n_yrs_diff[:, scen_rows(d)] .= n_yrs_above[:, scen_rows(d)] .- n_yrs_above[:, cf_idx(d)]
end

# Each performance map below is drawn once per location filter: filter 1 keeps only the reefs
# this scenario/DHW actually seeded, filter 2 keeps every reef seeded by at least one scenario
# (so all options share a common reef set, and reefs this option skipped are still shown).
filter_labels = ["seeded reefs", "all reefs"]
filter_names = ["seeded_reefs", "all_reefs"]
location_filters(d, opt_i) = [total_seeds[:, iv_col(d, opt_i)] .> 0, active_mask]
n_filters = length(filter_labels)

for (opt_i, scen_name) in enumerate(scenario_names), f in 1:n_filters
    # The seeded-reef filter depends on the DHW member, so resolve it per panel
    mask(d) = location_filters(d, opt_i)[f]

    fig_map = panel_map_figure(
        d -> (all_centroids[mask(d)], n_yrs_diff[mask(d), scen_row(d, opt_i)]),
        n_dhw;
        title_fn=d -> "Δ years above 20% coral cover — $scen_name vs counterfactual\n$(dhw_model_names[d]) — $(filter_labels[f])",
        cbar_label="Δ years above 20% coral cover (intervention − counterfactual)",
        colormap=Reverse(:vik25),
        colorrange=(-3, 3)
    )
    save(
        joinpath(
            pd_config["plot_output_path"],
            "n_yrs_above_target_$(scen_name)_$(filter_names[f]).png"
        ),
        fig_map;
        px_per_unit=2
    )
end

# Cumulative absolute cover difference map per scenario and location filter
# (one panel per DHW member)
cum_tac = dropdims(sum(m_tac; dims=1); dims=1)  # (n_locs, n_scens), km²·years
cum_tac_diff = Matrix{Float64}(undef, size(cum_tac))
for d in 1:n_dhw
    cum_tac_diff[:, scen_rows(d)] .= cum_tac[:, scen_rows(d)] .- cum_tac[:, cf_idx(d)]
end

# Symmetric colour range from the bulk of the deltas, so a few extreme reefs don't flatten the map
_cc_lim = quantile(abs.(vec(cum_tac_diff[active_mask, intervention_scens])), 0.98)
cum_tac_colorrange = iszero(_cc_lim) ? (-1.0, 1.0) :
    (-round(_cc_lim; sigdigits=2), round(_cc_lim; sigdigits=2))

for (opt_i, scen_name) in enumerate(scenario_names), f in 1:n_filters
    mask(d) = location_filters(d, opt_i)[f]

    fig_cc = panel_map_figure(
        d -> (all_centroids[mask(d)], cum_tac_diff[mask(d), scen_row(d, opt_i)]),
        n_dhw;
        title_fn=d -> "Δ cumulative coral cover — $scen_name vs counterfactual\n$(dhw_model_names[d]) — $(filter_labels[f])",
        cbar_label="Δ cumulative coral cover (km²·years, intervention − counterfactual)",
        colormap=Reverse(:vik),
        colorrange=cum_tac_colorrange
    )
    save(
        joinpath(
            pd_config["plot_output_path"],
            "cum_tac_diff_$(scen_name)_$(filter_names[f]).png"
        ),
        fig_cc;
        px_per_unit=2
    )
end

# Boxplots: performance metrics distribution per intervention scenario, dodged by DHW member
scen_labels = string.(scenario_names)
n_active = sum(active_mask)

x_scens = Int[]
dodge_dhw = Int[]
y_nyrs = Float64[]
y_cc = Float64[]
for opt_i in 1:length(scenario_names), d in 1:n_dhw
    s = scen_row(d, opt_i)
    append!(x_scens, fill(opt_i, n_active))
    append!(dodge_dhw, fill(d, n_active))
    append!(y_nyrs, n_yrs_diff[active_mask, s])
    append!(y_cc, cum_tac_diff[active_mask, s])
end
dhw_colors = Makie.wong_colors()[1:n_dhw]
box_colors = dhw_colors[dodge_dhw]

fig_box = Figure(; size=(1100, 450))
ax_nyrs = Axis(fig_box[1, 1];
    title="Δ years above 20% coral cover vs counterfactual",
    ylabel="Δ years above 20% coral cover",
    xticks=(1:length(scenario_names), scen_labels),
    xticklabelrotation=π / 4
)
boxplot!(ax_nyrs, x_scens, y_nyrs; dodge=dodge_dhw, color=box_colors)
hlines!(ax_nyrs, [0]; color=:black, linestyle=:dash, linewidth=1)

ax_cc = Axis(fig_box[1, 2];
    title="Δ cumulative coral cover vs counterfactual",
    ylabel="Δ cumulative coral cover (km²·years)",
    xticks=(1:length(scenario_names), scen_labels),
    xticklabelrotation=π / 4
)
boxplot!(ax_cc, x_scens, y_cc; dodge=dodge_dhw, color=box_colors)
hlines!(ax_cc, [0]; color=:black, linestyle=:dash, linewidth=1)

Legend(fig_box[1, 3],
    [PolyElement(; color=dhw_colors[d]) for d in 1:n_dhw],
    dhw_model_names;
    title="Global Climate Model", framevisible=false
)

save(
    joinpath(pd_config["plot_output_path"], "performance_metrics_boxplot.png"),
    fig_box;
    px_per_unit=2
)

# Options time-series plot
option_names = Symbol.(options.option_name)
all_names = vcat(option_names, [:unguided_intervention, :no_intervention])
intervention_names = all_names[1:(end - 1)]
# Groups index into a single DHW block (one figure per DHW member)
scen_groups = Dict{Symbol,BitVector}(
    name => BitVector((1:n_scens_per_dhw) .== i) for (i, name) in enumerate(all_names)
)
scen_groups_diff = Dict{Symbol,BitVector}(
    name => BitVector((1:n_options) .== i) for (i, name) in enumerate(intervention_names)
)

ts = string.(ADRIA.timesteps(rs))
tick_pos = collect(1:5:length(ts))
tick_lbl = ts[1:5:end]
(length(ts) - 1) % 5 != 0 &&
    (tick_pos = vcat(tick_pos, length(ts)); tick_lbl = vcat(tick_lbl, ts[end]))
xtick_vals = (tick_pos, tick_lbl)
xtick_rot = 2 / π

for (name, metric) in metrics, d in 1:n_dhw
    # This DHW block only; its counterfactual is the block's last column
    metric_dhw = ADRIA.DataCube(
        metric.data[:, scen_rows(d)];
        timesteps=ADRIA.timesteps(rs),
        scenarios=1:n_scens_per_dhw
    )
    metric_diff = ADRIA.DataCube(
        metric_dhw.data[:, 1:(end - 1)] .- metric_dhw.data[:, end];
        timesteps=ADRIA.timesteps(rs),
        scenarios=1:n_options
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

    save(
        joinpath(
            pd_config["plot_output_path"],
            "options_$(replace(lowercase(name), ' ' => '_'))_dhw$(dhw_scenarios[d]).png"
        ),
        f
    )
end
