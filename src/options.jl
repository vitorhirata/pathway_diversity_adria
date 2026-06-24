#=
Options experiment: compare five seeding strategy options across all DHW scenario members.

Design:
 5 options (from option_seed_preference) × N_dhw members = total scenarios
- Each option is held constant for the entire seeding period via option_ts,
  using pd_frequency = seed_years (one block = whole period)
- Only dhw_scenario varies; all other parameters are fixed
- SeedCriteriaWeights in the scenario spec are not relevant here —
  option_ts is what drives location selection in the guided path

To identify the option for scenario i after running:
    seed_year_start = Int(rs.inputs.seed_year_start[1])
    option = scens.option_ts[i][seed_year_start]  # e.g. :heat_stress
=#

include("src/common.jl")
using GeoMakie, GraphMakie, CairoMakie, NaturalEarth

RCP = "45"
seed_years = 30
dom = ADRIA.load_domain(pd_config["domain_path"], RCP;
    calib_params_fn=pd_config["coral_param_path"],
    # timeframe: seed_years + 2 (start seeding), 5 (extra years)
    timeframe=(2022, 2022 + seed_years + 2 + 5)
)
ms = ADRIA.model_spec(dom)

ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "FogCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "MCCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "Coral").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "GrowthAcceleration").fieldname)

N_seed_weights = (
    N_seed_TA=0.15, N_seed_CA=0.5, N_seed_CNA=0.0, N_seed_SM=0.35, N_seed_LM=0.0
)
N_seed_total = 1e7

ADRIA.fix_factor!(dom;
    # Environmental
    wave_scenario=1,
    dhw_scenario=11,
    # Intervention parameters
    guided=1,               # CoCoSo
    min_iv_locations=200,
    plan_horizon=5.0,
    # Alternative interventions
    N_mc_settlers=0,
    fogging=0.0,
    SRM=0.0,
    # Seeding parameters
    seed_year_start=2,
    seed_years=seed_years,
    seed_deployment_freq=1,
    seed_strategy=1,        # Periodic deployment
    seeding_devices_per_m2=5,
    a_adapt=5.0,
    a_adapt_ref=5,
    N_seed_TA=N_seed_total * N_seed_weights.N_seed_TA,
    N_seed_CA=N_seed_total * N_seed_weights.N_seed_CA,
    N_seed_CNA=N_seed_total * N_seed_weights.N_seed_CNA,
    N_seed_SM=N_seed_total * N_seed_weights.N_seed_SM,
    N_seed_LM=N_seed_total * N_seed_weights.N_seed_LM,
    # Depth
    depth_min=2.0,
    depth_offset=25.0
)

options = ADRIA.analysis.option_seed_preference(; include_weights=true)
scens = []
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
    push!(scens, scen)
end

# Unguided seeding
ADRIA.fix_factor!(dom; guided=0)
scen = ADRIA.sample(dom, 2)[1:1, :]
push!(scens, scen)

# No seeding
ADRIA.fix_factor!(dom;
    N_seed_TA=0,
    N_seed_CA=0,
    N_seed_CNA=0,
    N_seed_SM=0,
    N_seed_LM=0
)
scen = ADRIA.sample(dom, 2)[1:1, :]
push!(scens, scen)

scens = vcat(scens...)

rs = ADRIA.run_scenarios(dom, scens, RCP)

# Load scenario
# path = "Outputs/"
# rs = ADRIA.load_results(path)

# Reefs that received seeding in every intervention scenario
seed_start = Int(rs.inputs.seed_year_start[1])
n_seed_years = Int(rs.inputs.seed_years[1])
seed_ts = seed_start:(seed_start + n_seed_years - 1)
seed_per_reef_per_ts_scen = dropdims(
    sum(rs.seed_log[timesteps=seed_ts, scenarios=1:(nrow(scens) - 1)]; dims=:coral_id);
    dims=:coral_id
)

# always/never seeded: must hold across every timestep AND every scenario
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
    arrows2d!(ax, [arr_lon], [arr_lat0], [0.0], [arr_dlat];
        color=:black, tiplength=10.32, tipwidth=8.9)
    return text!(ax, arr_lon, arr_lat0 + arr_dlat + 0.12;
        text="N", align=(:center, :bottom), fontsize=13, font=:bold)
end

#=
# Selected locations GIF per intervention scenario
plottable_gif = GeoMakie.to_multipoly(dom.loc_data[:, :geometry])
for (scen_idx, scen_name) in enumerate(scenario_names)
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
        joinpath(pd_config["plot_output_path"], "seeding_map_$(scen_name).gif"),
        eachindex(ts_labels);
        framerate=3
    ) do i
        seeded_points[] = all_centroids[seed_per_reef_per_ts_scen[timesteps=i, scenarios=scen_idx] .> 0]
        title_obs[] = "$scen_name — Year: $(ts_labels[i])"
    end
end
=#

# Seeding frequency histograms per intervention scenario
fig_hist = Figure(; size=(length(scenario_names) * 350, 400))
for (scen_idx, scen_name) in enumerate(scenario_names)
    seeded_binary = seed_per_reef_per_ts_scen[scenarios=scen_idx] .> 0
    seeding_freq = vec(mean(seeded_binary; dims=1)) .* 100
    ax = Axis(
        fig_hist[1, scen_idx];
        xlabel="Seeding frequency (%)",
        ylabel=scen_idx == 1 ? "Number of reefs" : "",
        title=string(scen_name),
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
    title="Seeding frequency — all intervention scenarios",
    xticks=0:10:100
)
hist!(ax_hist_all, seeding_freq_all; bins=0:5:100)
save(joinpath(pd_config["plot_output_path"], "seeding_frequency_all.png"), fig_hist_all)

cf_idx = nrow(scens)  # last scenario = no_intervention
active_mask = .!never_seeded

# Total seeds per location map per intervention scenario
total_seeds = Array(
    dropdims(
        sum(rs.seed_log[scenarios=1:(nrow(scens) - 1)]; dims=(:timesteps, :coral_id));
        dims=(:timesteps, :coral_id)
    )
)  # (n_locs, n_intervention_scens)
#seed_colorrange = (0.0, Float64(maximum(total_seeds)))
seed_colorrange = (0.0, 1.2e7)

for (scen_idx, scen_name) in enumerate(scenario_names)
    scen_no_seeds = total_seeds[:, scen_idx] .== 0
    scen_has_seeds = .!scen_no_seeds

    fig_seeds = Figure()
    ga_seeds = GeoAxis(
        fig_seeds[1, 1];
        dest="+proj=longlat +datum=WGS84",
        limits=(_gbr_lon_min, _gbr_lon_max, _gbr_lat_min, _gbr_lat_max),
        title="Total deployed coral per location — $scen_name",
        xgridcolor=(:gray, 0.15),
        ygridcolor=(:gray, 0.15),
        xticklabelsize=7,
        yticklabelsize=7
    )
    poly!(
        ga_seeds,
        _ne_land.geometry;
        color=RGBf(0.93, 0.91, 0.87),
        strokewidth=0.5,
        strokecolor=:gray40
    )
    scatter!(ga_seeds, all_centroids[scen_no_seeds]; color=:gray80, markersize=4, alpha=0.5)
    seeds_vals = total_seeds[scen_has_seeds, scen_idx]
    order = sortperm(seeds_vals)
    scatter!(ga_seeds, all_centroids[scen_has_seeds][order]; color=seeds_vals[order],
        colormap=:viridis, colorrange=seed_colorrange, markersize=4, alpha=0.7)
    Colorbar(fig_seeds[1, 2];
        colorrange=seed_colorrange,
        colormap=:viridis,
        label="Total deployed coral",
        height=Relative(0.65)
    )
    _gbr_annotations!(ga_seeds)
    rowsize!(fig_seeds.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig_seeds)
    save(
        joinpath(pd_config["plot_output_path"], "total_seeds_$(scen_name).png"),
        fig_seeds;
        px_per_unit=2
    )
end

# n_yrs_above_target scatter map per scenario
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

for (scen_idx, scen_name) in enumerate(scenario_names)
    scen_no_seeds = total_seeds[:, scen_idx] .== 0
    scen_has_seeds = .!scen_no_seeds

    location_filter = active_mask # scen_has_seeds or active_mask
    diff = n_yrs_above[location_filter, scen_idx] .- n_yrs_above[location_filter, cf_idx]
    order = sortperm(diff)

    fig_map = Figure()
    ga_map = GeoAxis(
        fig_map[1, 1];
        dest="+proj=longlat +datum=WGS84",
        limits=(_gbr_lon_min, _gbr_lon_max, _gbr_lat_min, _gbr_lat_max),
        title="Δ years above 20% coral cover — $scen_name vs couterfactual",
        xgridcolor=(:gray, 0.15),
        ygridcolor=(:gray, 0.15),
        xticklabelsize=7,
        yticklabelsize=7
    )
    diff_all = n_yrs_above[active_mask, scen_idx] .- n_yrs_above[active_mask, cf_idx]
    ann_med = round(median(diff_all); digits=1)
    ann_p25 = round(quantile(diff_all, 0.25); digits=1)
    ann_p75 = round(quantile(diff_all, 0.75); digits=1)
    ann_lat = _gbr_lat_max + 0.25 * (_gbr_lat_min - _gbr_lat_max)

    poly!(
        ga_map,
        _ne_land.geometry;
        color=RGBf(0.93, 0.91, 0.87),
        strokewidth=0.5,
        strokecolor=:gray40
    )
    #scatter!(ga_map, all_centroids[scen_no_seeds]; color=:gray80, markersize=3, alpha=0.5)
    scatter!(ga_map, all_centroids[location_filter][order]; color=diff[order],
        colormap=Reverse(:vik25),
        colorrange=(-10, 10), markersize=4, alpha=0.7)
    Colorbar(fig_map[1, 2];
        colorrange=(-10, 10),
        colormap=Reverse(:vik25),
        label="Δ years above 20% coral cover (intervention − counterfactual)",
        height=Relative(0.65)
    )
    _gbr_annotations!(ga_map)
    text!(ga_map, _gbr_lon_max - 0.2, ann_lat;
        text="median=$(ann_med) ($(ann_p25)–$(ann_p75))",
        align=(:right, :center), fontsize=8, color=:black)
    rowsize!(fig_map.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig_map)
    save(
        joinpath(pd_config["plot_output_path"], "n_yrs_above_target_$(scen_name).png"),
        fig_map;
        px_per_unit=2
    )
end

# Cumulative functional diversity difference map per scenario
fd_data = Array(ADRIA.metrics.coral_evenness(rs))
cf_fd = fd_data[:, :, cf_idx]
cum_fd_diff = dropdims(sum(fd_data .- cf_fd; dims=1); dims=1)  # (n_locs, n_scens)

for (scen_idx, scen_name) in enumerate(scenario_names)
    scen_no_seeds = total_seeds[:, scen_idx] .== 0
    scen_has_seeds = .!scen_no_seeds

    location_filter = active_mask # scen_has_seeds or active_mask
    diff_fd = cum_fd_diff[location_filter, scen_idx]
    order_fd = sortperm(diff_fd)

    fig_fd = Figure()
    ga_fd = GeoAxis(
        fig_fd[1, 1];
        dest="+proj=longlat +datum=WGS84",
        limits=(_gbr_lon_min, _gbr_lon_max, _gbr_lat_min, _gbr_lat_max),
        title="Δ cumulative coral evenness — $scen_name vs counterfactual",
        xgridcolor=(:gray, 0.15),
        ygridcolor=(:gray, 0.15),
        xticklabelsize=7,
        yticklabelsize=7
    )
    diff_fd_all = cum_fd_diff[active_mask, scen_idx]
    ann_med_fd = round(median(diff_fd_all); digits=2)
    ann_p25_fd = round(quantile(diff_fd_all, 0.25); digits=2)
    ann_p75_fd = round(quantile(diff_fd_all, 0.75); digits=2)
    ann_lat_fd = _gbr_lat_max + 0.25 * (_gbr_lat_min - _gbr_lat_max)

    poly!(
        ga_fd,
        _ne_land.geometry;
        color=RGBf(0.93, 0.91, 0.87),
        strokewidth=0.5,
        strokecolor=:gray40
    )
    #scatter!(ga_fd, all_centroids[scen_no_seeds]; color=:gray80, markersize=3, alpha=0.5)
    scatter!(ga_fd, all_centroids[location_filter][order_fd]; color=diff_fd[order_fd],
        colormap=Reverse(:vik), colorrange=(-10, 10), markersize=4, alpha=0.7)
    Colorbar(fig_fd[1, 2];
        colorrange=(-10, 10),
        colormap=Reverse(:vik),
        label="Δ cumulative coral evenness (intervention − counterfactual)",
        height=Relative(0.65)
    )
    _gbr_annotations!(ga_fd)
    text!(ga_fd, _gbr_lon_max - 0.2, ann_lat_fd;
        text="median=$(ann_med_fd) ($(ann_p25_fd)–$(ann_p75_fd))",
        align=(:right, :center), fontsize=8, color=:black)
    rowsize!(fig_fd.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig_fd)
    save(
        joinpath(pd_config["plot_output_path"], "cum_fd_diff_$(scen_name).png"),
        fig_fd;
        px_per_unit=2
    )
end

# Boxplots: performance metrics distribution per intervention scenario
scen_labels = string.(scenario_names)
n_active = sum(active_mask)
x_scens = vcat([fill(i, n_active) for i in 1:length(scenario_names)]...)

y_nyrs = vcat(
    [
        n_yrs_above[active_mask, s] .- n_yrs_above[active_mask, cf_idx]
        for s in 1:length(scenario_names)
    ]...
)
y_fd = vcat([cum_fd_diff[active_mask, s] for s in 1:length(scenario_names)]...)

fig_box = Figure(; size=(1000, 450))
ax_nyrs = Axis(fig_box[1, 1];
    title="Δ years above 20% coral cover vs counterfactual",
    ylabel="Δ years above 20% coral cover",
    xticks=(1:length(scenario_names), scen_labels),
    xticklabelrotation=π / 4
)
boxplot!(ax_nyrs, x_scens, Float64.(y_nyrs))
hlines!(ax_nyrs, [0]; color=:black, linestyle=:dash, linewidth=1)

ax_fd = Axis(fig_box[1, 2];
    title="Δ cumulative coral evenness vs counterfactual",
    ylabel="Δ cumulative coral evenness",
    xticks=(1:length(scenario_names), scen_labels),
    xticklabelrotation=π / 4
)
boxplot!(ax_fd, x_scens, y_fd)
hlines!(ax_fd, [0]; color=:black, linestyle=:dash, linewidth=1)

save(
    joinpath(pd_config["plot_output_path"], "performance_metrics_boxplot.png"),
    fig_box;
    px_per_unit=2
)

# Options time-series plot
option_names = Symbol.(options.option_name)
all_names = vcat(option_names, [:unguided_intervention, :no_intervention])
intervention_names = all_names[1:(end - 1)]
scen_groups = Dict{Symbol,BitVector}(
    name => BitVector((1:nrow(scens)) .== i) for (i, name) in enumerate(all_names)
)
scen_groups_diff = Dict{Symbol,BitVector}(
    name => BitVector((1:(nrow(scens) - 1)) .== i) for
    (i, name) in enumerate(intervention_names)
)

ts = string.(ADRIA.timesteps(rs))
tick_pos = collect(1:5:length(ts))
tick_lbl = ts[1:5:end]
(length(ts) - 1) % 5 != 0 &&
    (tick_pos = vcat(tick_pos, length(ts)); tick_lbl = vcat(tick_lbl, ts[end]))
xtick_vals = (tick_pos, tick_lbl)
xtick_rot = 2 / π

for (name, metric) in metrics
    metric_diff = ADRIA.DataCube(
        metric.data[:, 1:(end - 1)] .- metric.data[:, end];
        timesteps=ADRIA.timesteps(rs),
        scenarios=1:(nrow(scens) - 1)
    )

    f = Figure(; size=(3200, 800))
    g1 = f[1, 1] = GridLayout()
    g2 = f[1, 2] = GridLayout()

    ax1 = Axis(g1[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot, title=name)
    ADRIA.viz.scenarios!(g1, ax1, metric, scen_groups;
        opts=Dict{Symbol,Any}(
            :legend_labels => all_names, :legend => false, :histogram => false
        ))
    ADRIA.viz.scenarios_legend!(g1[1, 0], scen_groups, metric;
        opts=Dict{Symbol,Any}(:legend_labels => all_names),
        legend_opts=Dict{Symbol,Any}(:padding => (4, 4, 4, 4))
    )

    ax2 = Axis(g2[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot,
        title="$name - counterfactual", ylabel="$name - counterfactual")
    ADRIA.viz.scenarios!(g2, ax2, metric_diff, scen_groups_diff;
        opts=Dict{Symbol,Any}(
            :legend_labels => intervention_names, :legend => false, :histogram => false
        ))

    save(
        joinpath(
            pd_config["plot_output_path"],
            "options_$(replace(lowercase(name), ' ' => '_')).png"
        ),
        f
    )
end
