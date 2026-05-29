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
using GeoMakie, GraphMakie, WGLMakie

RCP = "45"
dom = ADRIA.load_domain(pd_config["domain_path"], RCP; calib_params_fn=pd_config["coral_param_path"])
ms = ADRIA.model_spec(dom)

ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "FogCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "MCCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "Coral").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "GrowthAcceleration").fieldname)

N_seed_weights = (N_seed_TA=0.15, N_seed_CA=0.5, N_seed_CNA=0.0, N_seed_SM=0.35, N_seed_LM=0.0)
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
    seed_years=30,
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

options = ADRIA.analysis.option_seed_preference(include_weights=true)
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
    N_seed_LM=0,
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
    sum(rs.seed_log[timesteps=seed_ts, scenarios=1:(nrow(scens) - 1)]; dims=:coral_id),
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
s_juves = ADRIA.metrics.scenario_relative_juveniles(_aj[locations=selected_locations].data, _k_area)
metrics = Dict(
    "Total absolute cover"    => s_tac,
    "Relative Shelter Volume" => s_rsv,
    "Relative Juveniles"      => s_juves,
    "Coral Evenness"          => s_even
)

# Selected locations GIF per intervention scenario
ts_labels = ADRIA.timesteps(rs)[seed_ts]
all_centroids = ADRIA.centroids(dom.loc_data)
plottable_gif = GeoMakie.to_multipoly(dom.loc_data[:, :geometry])
scenario_names = vcat(options.option_name, [:unguided])

for (scen_idx, scen_name) in enumerate(scenario_names)
    seeded_points = Observable(Point2f[])
    title_obs = Observable("$scen_name — Year: $(ts_labels[1])")

    fig_gif = Figure(; size=(650, 920), figure_padding=5)
    ga_gif = GeoAxis(
        fig_gif[1, 1];
        dest="+proj=latlong +datum=WGS84",
        title=title_obs,
        titlesize=20,
        aspect=DataAspect(),
        limits=((141.8, 153.7), (-25.2, -9.8)),
        xgridwidth=0.5,
        ygridwidth=0.5,
    )
    poly!(ga_gif, plottable_gif; color=:gray80)
    scatter!(ga_gif, seeded_points; color=:red, markersize=4)

    record(fig_gif, joinpath(pd_config["plot_output_path"], "seeding_map_$(scen_name).gif"), eachindex(ts_labels); framerate=3) do i
        seeded_points[] = all_centroids[seed_per_reef_per_ts_scen[timesteps=i, scenarios=scen_idx] .> 0]
        title_obs[] = "$scen_name — Year: $(ts_labels[i])"
    end
end

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
        xticks=0:20:100,
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
    xticks=0:10:100,
)
hist!(ax_hist_all, seeding_freq_all; bins=0:5:100)
save(joinpath(pd_config["plot_output_path"], "seeding_frequency_all.png"), fig_hist_all)

# Options time-series plot
option_names = Symbol.(options.option_name)
all_names = vcat(option_names, [:unguided_intervention, :no_intervention])
intervention_names = all_names[1:(end-1)]
scen_groups = Dict{Symbol,BitVector}(
    name => BitVector((1:nrow(scens)) .== i) for (i, name) in enumerate(all_names)
)
scen_groups_diff = Dict{Symbol,BitVector}(
    name => BitVector((1:(nrow(scens) - 1)) .== i) for (i, name) in enumerate(intervention_names)
)

ts = string.(ADRIA.timesteps(rs))
tick_pos = collect(1:5:length(ts))
tick_lbl = ts[1:5:end]
(length(ts) - 1) % 5 != 0 && (tick_pos = vcat(tick_pos, length(ts)); tick_lbl = vcat(tick_lbl, ts[end]))
xtick_vals = (tick_pos, tick_lbl)
xtick_rot = 2 / π

for (name, metric) in metrics
    metric_diff = ADRIA.DataCube(
        metric.data[:, 1:(end-1)] .- metric.data[:, end];
        timesteps=ADRIA.timesteps(rs),
        scenarios=1:(nrow(scens) - 1)
    )

    f = Figure(; size=(3200, 800))
    g1 = f[1, 1] = GridLayout()
    g2 = f[1, 2] = GridLayout()

    ax1 = Axis(g1[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot, title=name)
    ADRIA.viz.scenarios!(g1, ax1, metric, scen_groups;
        opts=Dict{Symbol,Any}(:legend_labels => all_names, :legend => false, :histogram => false))
    ADRIA.viz.scenarios_legend!(g1[1, 0], scen_groups, metric;
        opts=Dict{Symbol,Any}(:legend_labels => all_names),
        legend_opts=Dict{Symbol,Any}(:padding => (4, 4, 4, 4))
    )

    ax2 = Axis(g2[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot,
        title="$name - counterfactual", ylabel="$name - counterfactual")
    ADRIA.viz.scenarios!(g2, ax2, metric_diff, scen_groups_diff;
        opts=Dict{Symbol,Any}(:legend_labels => intervention_names, :legend => false, :histogram => false))

    save(joinpath(pd_config["plot_output_path"], "options_$(replace(lowercase(name), ' ' => '_')).png"), f)
end
