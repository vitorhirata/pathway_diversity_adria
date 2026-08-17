include("src/common.jl")
using GeoMakie, GraphMakie, CairoMakie

rcps = ["26", "45", "70"]
scen_ids = 1:11
start_year = 2022
last_year = 2100
dom = ADRIA.load_domain(pd_config["domain_path"], rcps[1])

for rcp in rcps
    dom = ADRIA.switch_RCPs!(dom, rcp)

    fig = Figure(size=(900, 700))
    ax = Axis(fig[1, 1];
        title="DHW distribution — RCP $(rcp) Years $(start_year)-$(last_year)",
        xlabel="Degree Heating Weeks (°C-weeks)",
        ylabel="DHW scenario (GCM model)",
        yticks=(collect(scen_ids), dom.dhw_scens.properties["model_names"]),
        limits=(0, 70, nothing, nothing)
    )

    for scen_id in scen_ids
        dhw_scens = dom.dhw_scens[
            scenarios=scen_id, mcb_durations=ADRIA.At(0), albedo=1,
            timesteps=ADRIA.At(start_year:last_year)
        ]

        # flatten across timesteps and locations into a single vector
        y = vec(Array(dhw_scens))
        boxplot!(ax, fill(scen_id, length(y)), y; orientation=:horizontal)
    end

    save(
        joinpath(pd_config["plot_output_path"], "dhw_distribution_rcp$(rcp)_years_$(start_year)-$(last_year).png"),
        fig; px_per_unit=2
    )
    @info "Saved dhw_distribution_rcp$(rcp)_years_$(start_year)-$(last_year).png"
end

# ── DHW time-series: median over locations per (GCM, RCP), decade P10-P90 whiskers ──

dhw_scenarios = [5, 11]
years = collect(start_year:last_year)
decade_years = start_year:10:last_year

# centered rolling average (window shrinks at the edges to keep the series length)
function rolling_mean(x, w)
    n = length(x)
    half = w ÷ 2
    return [mean(x[max(1, i - half):min(n, i + half)]) for i in 1:n]
end

# distinct color per (GCM, RCP) combination
n_combos = length(dhw_scenarios) * length(rcps)
combo_colors = cgrad(:tab10, n_combos; categorical=true)

fig = Figure(size=(1000, 700))
ax = Axis(fig[1, 1];
    title="DHW median over locations — Years $(start_year)-$(last_year)",
    xlabel="Year",
    ylabel="Degree Heating Weeks (°C-weeks)"
)

c = 0
for rcp in rcps[2:3]
    dom = ADRIA.switch_RCPs!(dom, rcp)
    model_names = dom.dhw_scens.properties["model_names"]

    for scen_id in dhw_scenarios
        global c += 1
        color = combo_colors[c]

        dhw_scens = dom.dhw_scens[
            scenarios=scen_id, mcb_durations=ADRIA.At(0), albedo=1,
            timesteps=ADRIA.At(start_year:last_year)
        ]
        arr = Array(dhw_scens)  # (timesteps, locations)

        # median over locations per timestep, smoothed with a 5-year rolling average
        med = rolling_mean(vec(median(arr; dims=2)), 5)
        lines!(ax, years, med;
            color=color, label="RCP $(rcp) — $(model_names[scen_id])"
        )

        #=
        # decade markers: median + P10-P90 whiskers over locations
        dec_idx = [findfirst(==(yr), years) for yr in decade_years]
        dec_med = [median(arr[i, :]) for i in dec_idx]
        dec_p10 = [quantile(arr[i, :], 0.10) for i in dec_idx]
        dec_p90 = [quantile(arr[i, :], 0.90) for i in dec_idx]
        errorbars!(ax, collect(decade_years), dec_med, dec_med .- dec_p10, dec_p90 .- dec_med;
            color=color, whiskerwidth=8
        )
        scatter!(ax, collect(decade_years), dec_med; color=color, markersize=8)
        =#
    end
end

axislegend(ax; position=:lt, framevisible=false)

save(
    joinpath(pd_config["plot_output_path"], "dhw_timeseries_years_$(start_year)-$(last_year).png"),
    fig; px_per_unit=2
)
@info "Saved dhw_timeseries_years_$(start_year)-$(last_year).png"
