#=
Reef polygon overview map.

Loads the ADRIA domain and draws the polygon of every reef in `dom.loc_data` (the `:geometry`
column) on a Great Barrier Reef map, in a single uniform colour. A visual/orientation aid
showing the shape and location of all reefs in the data package.

Reuses the shared GBR map style from `visualization/static_options.jl` (`_ne_land` land fill,
`_gbr_annotations!` city labels + scale bar + north arrow, and the `_gbr_*` extent constants).
Only those pieces are used, so the label/scenario globals that file expects at call time are
not needed here. Run from the repo root:

    julia --project src/one-off_scripts/reef_polygons_map.jl
=#

include("src/common.jl")
using GeoMakie, CairoMakie, NaturalEarth
using Statistics
include("src/visualization/robustness.jl")      # defines `_sci`, required before the next include
include("src/visualization/static_options.jl")  # reuse `_ne_land`, `_gbr_annotations!`, `_gbr_*` extents

# ── Domain ────────────────────────────────────────────────────────────────────
# Minimal load: geometry is independent of RCP, calibration, and factor settings.
dom = ADRIA.load_domain(pd_config["domain_path"], "45")
@info "Loaded domain" n_reefs = nrow(dom.loc_data) columns = names(dom.loc_data)

# ── Map ─────────────────────────────────────────────────────────────────────
reef_polys = GeoMakie.to_multipoly(dom.loc_data[:, :geometry])

fig = Figure(; size=(1200, 900))

# 2° ticks + graticule, labelled every other degree. A plain Axis is used instead of GeoAxis
# (the map is equirectangular, so polygons render identically in lon/lat) because GeoAxis
# auto-filters ticks near the frame — dropping most labels for a sparse 2° set. Same approach
# as `plot_total_seeds` in `visualization/static_options.jl`.
grid_lons = ceil(Int, _gbr_lon_min):floor(Int, _gbr_lon_max)
grid_lats = ceil(Int, _gbr_lat_min):floor(Int, _gbr_lat_max)
label_lons = collect(filter(iseven, grid_lons))
label_lats = collect(filter(iseven, grid_lats))

ax = Axis(
    fig[1, 1];
    limits=(_gbr_lon_min, _gbr_lon_max, _gbr_lat_min, _gbr_lat_max),
    title="",
    xgridvisible=false,
    ygridvisible=false,
    xticks=(label_lons, ["$(l)°E" for l in label_lons]),
    yticks=(label_lats, ["$(abs(l))°S" for l in label_lats]),
    xticklabelsize=14,
    yticklabelsize=14
)
hidespines!(ax)  # drop the axis frame box

# Land fill (same style as `panel_map_figure`)
poly!(
    ax,
    _ne_land.geometry;
    color=RGBf(0.93, 0.91, 0.87),
    strokewidth=0.5,
    strokecolor=:gray40
)

# 2° graticule over the whole extent
for lon in label_lons
    lines!(ax, [lon, lon], [_gbr_lat_min, _gbr_lat_max]; color=(:gray, 0.15), linewidth=0.5)
end
for lat in label_lats
    lines!(ax, [_gbr_lon_min, _gbr_lon_max], [lat, lat]; color=(:gray, 0.15), linewidth=0.5)
end

# All reef polygons — single uniform colour
poly!(
    ax,
    reef_polys;
    color=(:teal, 0.75),
    strokecolor=(:black, 0.3),
    strokewidth=0.2
)

_gbr_annotations!(
    ax;
    label_min_lon=145.0, label_min_lat=-20.0, scalerank=4,
    city_fontsize=14,                          # larger town labels
    scale_lon0=7.0, scale_km=200.0, scale_fontsize = 12, # scale bar: 200 km (double the default)
    arrow_lon=_gbr_lon_min + 3.5, arrow_lat=_gbr_lat_max - 2.0  # tweak to reposition the N arrow
)

# ── Australia locator inset (inside the main map, bottom-left) ────────────────
# Australia-only map (other countries/islands clipped out) with a red box marking the GBR
# extent shown in the main panel. Positioned inside the main map with its left edge after 142°E.
_aus_lon_min, _aus_lon_max = 112.0, 154.0
_aus_lat_min, _aus_lat_max = -44.0, -8.0  # top margin so the GBR box (to -9.8°) isn't clipped

# Australia country polygon only
_ne_countries = naturalearth("admin_0_countries", 10)
_aus_geom = [f.geometry for f in _ne_countries if get(f.properties, :ADMIN, "") == "Australia"]

inset_ax = GeoAxis(
    fig[1, 1];
    width=Relative(0.35),
    height=Relative(0.34),
    halign=0.1,       # fraction of free space → left edge sits just east of 142°E
    valign=0.0,      # nudged up off the bottom so the x-axis labels stay clear
    dest="+proj=longlat +datum=WGS84",
    limits=(_aus_lon_min, _aus_lon_max, _aus_lat_min, _aus_lat_max)
)
hidedecorations!(inset_ax)

# White background to match the main map's ocean (GeoAxis has no backgroundcolor attribute)
poly!(
    inset_ax,
    Rect2f(_aus_lon_min, _aus_lat_min, _aus_lon_max - _aus_lon_min, _aus_lat_max - _aus_lat_min);
    color=RGBf(0.93, 0.91, 0.87)
)

poly!(
    inset_ax,
    GeoMakie.to_multipoly(_aus_geom);
    color=RGBf(0.97, 0.71, 0.2),
    strokewidth=0.3,
    strokecolor=:gray40
)

# Red box outlining the GBR extent (no fill)
lines!(
    inset_ax,
    [_gbr_lon_min, _gbr_lon_max, _gbr_lon_max, _gbr_lon_min, _gbr_lon_min],
    [_gbr_lat_min, _gbr_lat_min, _gbr_lat_max, _gbr_lat_max, _gbr_lat_min];
    color=:red,
    linewidth=1.5
)

# ── DHW time-series panel (right) ─────────────────────────────────────────────
# Mean-over-reefs DHW trajectory with a 95% percentile band across reefs, for three DHW
# scenario members, labelled by climate-model name, against a 2007–2015 historical baseline.
dhw_scenarios = [2, 7, 10]
dhw_model_names = dom.dhw_scens.properties["model_names"][dhw_scenarios]
years = collect(dom.env_layer_md.timeframe)
ci_level = 0.95
α = (1 - ci_level) / 2
plot_start_year = 2022  # left x-limit of the DHW panel
plot_end_year = 2050    # right x-limit of the DHW panel
xtick_years = (cld(plot_start_year, 5) * 5):5:plot_end_year  # ticks every 5 years

dhw_ax = Axis(
    fig[1, 2];
    xlabel="Year",
    ylabel="DHW (°C-week)",
    ylabelsize=16,
    xlabelsize=16,
    xticklabelsize=15,
    yticklabelsize=15,
    xticks=xtick_years,
    limits=(plot_start_year, plot_end_year, 0, 40)
)
hidespines!(dhw_ax, :t, :r)

colors = Makie.wong_colors()[1:length(dhw_scenarios)]

# This is an MCB data package, so `dhw_scens` is 5-D (timesteps, locations, scenarios,
# mcb_durations, albedo). Take the no-MCB baseline slice (duration index 1, lowest albedo),
# matching ADRIA's `dhw_baseline` in `scenario.jl`, to get a plain (timesteps, locations) array.
dhw_baseline(s) = dom.dhw_scens[
    scenarios=s, mcb_durations=1, albedo=At(dom.dhw_scens.albedo[1])
]

plot_handles = []  # one mean-line handle per scenario, for the legend
for (i, s) in enumerate(dhw_scenarios)
    data_s = dhw_baseline(s)  # (timesteps, locations)
    # Iterate a plain range (not `axes`, which is dimensional here) so the comprehensions
    # return plain `Vector{Float64}` that Makie can plot directly.
    means = Float64[mean(data_s[t, :]) for t in 1:size(data_s, 1)]
    lower = Float64[quantile(vec(data_s[t, :]), α) for t in 1:size(data_s, 1)]
    upper = Float64[quantile(vec(data_s[t, :]), 1 - α) for t in 1:size(data_s, 1)]

    band!(dhw_ax, years, lower, upper; color=(colors[i], 0.25))
    ln = lines!(dhw_ax, years, means; color=colors[i], linewidth=2)
    push!(plot_handles, ln)
end

# Historical baseline: mean DHW over 2007–2015 across all reefs and scenarios (identical for
# every member), drawn as a dashed horizontal line.
hist_idx = findall(y -> 2007 <= y <= 2015, years)
hist_baseline = mean([mean(dhw_baseline(s)[hist_idx, :]) for s in dhw_scenarios])
hb = hlines!(dhw_ax, hist_baseline; color=:black, linestyle=:dash, linewidth=1.5)

# Summary: mean DHW over the plotted period, across all reefs, per scenario — both raw and as
# a percent change relative to the historical baseline.
plot_idx = findall(y -> plot_start_year <= y <= plot_end_year, years)
mean_dhw_period = [mean(dhw_baseline(s)[plot_idx, :]) for s in dhw_scenarios]
pct_vs_hist = (mean_dhw_period ./ hist_baseline .- 1) .* 100
@info "Mean DHW over plotted period ($(plot_start_year)–$(plot_end_year)), all reefs" historical_baseline = round(hist_baseline; digits=2) model = dhw_model_names mean_dhw = round.(mean_dhw_period; digits=2) pct_change_vs_baseline = round.(pct_vs_hist; digits=1)

axislegend(
    dhw_ax,
    vcat(plot_handles, hb),
    vcat(string.(dhw_model_names), "Historical baseline (2007–2015)");
    position=:lt,
    labelsize=16,
    framevisible=false
)

# ── Panel labels ──────────────────────────────────────────────────────────────
for (col, lab) in ((1, "(a)"), (2, "(b)"))
    Label(
        fig[1, col, TopLeft()], lab;
        font=:bold, fontsize=18, halign=:left, valign=:bottom, padding=(0, 0, 5, 0)
    )
end

# ── Layout sizing ─────────────────────────────────────────────────────────────
# Cos-corrected aspect so the map isn't horizontally stretched (plain Axis has no projection).
map_aspect = (_gbr_lon_max - _gbr_lon_min) * cosd((_gbr_lat_min + _gbr_lat_max) / 2) /
             (_gbr_lat_max - _gbr_lat_min)
rowsize!(fig.layout, 1, Fixed(820))
colsize!(fig.layout, 1, Aspect(1, map_aspect))
colsize!(fig.layout, 2, Fixed(520))  # DHW panel width (col 1 is pinned to the map aspect)
resize_to_layout!(fig)

# ── Save ────────────────────────────────────────────────────────────────────
save(joinpath(pd_config["plot_output_path"], "reef_polygons_dhw.png"), fig; px_per_unit=2)
@info "Saved reef_polygons_dhw.png"
