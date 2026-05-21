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

using Revise
using Infiltrator
using Statistics, DataFrames
using ADRIA

RME_path = "/home/vitor/Code/ADRIA.jl/DataPackages/GBR_MCB_GBR_2026-03-30_v080/"
calib_path = "/home/vitor/Code/ADRIA.jl/DataPackages/calibrated_params.nc"
RCP = "45"

dom = ADRIA.load_domain(RME_path, RCP; calib_params_fn=calib_path)
ms = ADRIA.model_spec(dom)

ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "FogCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "MCCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "Coral").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "GrowthAcceleration").fieldname)

ADRIA.fix_factor!(dom;
    # Environmental
    wave_scenario=1,
    dhw_scenario=11,
    # Intervention
    guided=1,               # CoCoSo
    N_mc_settlers=0,
    seeding_devices_per_m2=5,
    fogging=0.0,
    SRM=0.0,
    a_adapt=0.0,
    a_adapt_ref=0.0,
    seed_years=30,
    seed_deployment_freq=1,
    plan_horizon=20.0,
    seed_year_start=2,
    seed_strategy=1,        # Periodic deployment
    min_iv_locations=1000,
    # 1e8 corals total split equally across 5 species
    N_seed_TA=2e7,
    N_seed_CA=2e7,
    N_seed_CNA=2e7,
    N_seed_SM=2e7,
    N_seed_LM=2e7,
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

# Reefs that received seeding in every intervention scenario (across timesteps 2-32)
seed_per_reef_per_scen = dropdims(
    sum(rs.seed_log[timesteps=2:32, scenarios=1:(n_scens - 1)]; dims=(:timesteps, :coral_id)),
    dims=(:timesteps, :coral_id)
)
always_seeded = vec(all(seed_per_reef_per_scen.data .> 0; dims=2))
never_seeded = vec(all(seed_per_reef_per_scen.data .== 0; dims=2))

selected_locations = .!(always_seeded .| never_seeded) # Remove reefs that always or never have seeding
selected_locations = findall(selected_locations)

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

using GeoMakie, GraphMakie, WGLMakie

# Grouping setup
option_names = Symbol.(options.option_name)
all_names = vcat(option_names, [:unguided_intervention, :no_intervention])
intervention_names = all_names[1:(end-1)]
n_scens = nrow(scens)
scen_groups = Dict{Symbol,BitVector}(
    name => BitVector((1:n_scens) .== i) for (i, name) in enumerate(all_names)
)
scen_groups_diff = Dict{Symbol,BitVector}(
    name => BitVector((1:(n_scens - 1)) .== i) for (i, name) in enumerate(intervention_names)
)

# Shared x-axis ticks
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
        scenarios=1:(n_scens - 1)
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

    save("options_$(replace(lowercase(name), ' ' => '_')).png", f)
end
