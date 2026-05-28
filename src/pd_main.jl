#=
Pathway diversity main simulation.
Main parameters varied:
- RCP
- dhw_scenario
- min_iv_locations
- Number of seeds
=#

include("src/common.jl")
using WGLMakie, GeoMakie, GraphMakie

RCP = "26" # RCP 26, 45, 70
seed_years = 20
dom = ADRIA.load_domain(
    pd_config["domain_path"], RCP;
    calib_params_fn=pd_config["coral_param_path"], timeframe=(2022, 2022 + 20 + seed_years + 2)
)
ms = ADRIA.model_spec(dom)

ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "FogCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "MCCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "Coral").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "GrowthAcceleration").fieldname)

# Generate scenarios
ADRIA.fix_factor!(dom;
    #Seeding params
    seed_year_start=2,     # Start as soon as possible. Possible parameter to vary
    seed_years=seed_years, # Based on pathway diversity analysis time
    seed_deployment_freq=1,# Lower bound of distribution. Seed every year.
    seeding_devices_per_m2=5,
    seed_strategy=1,        # Periodic deployment
    a_adapt=5.0,
    a_adapt_ref=5,
    # Interventions params
    plan_horizon=5.0,     # Upper bound of distribution
    guided=1,              # CoCoSo. Better performance based on initial analysis
    #Decision params. Chosen to allow more locations receive intervention
    depth_min=2.0,         # Lower bound of distribution
    depth_offset=25.0,     # Upper bound of distribution
    #Other interventions (no fogging, shading or moving corals)
    fogging=0.0,
    SRM=0.0,
    N_mc_settlers=0,
    #Environmental params
    wave_scenario=1
)

dhw_scenarios = [5, 7, 11]
n_seed_locations = [
    [1e2, 15],
    [1e6, 200],
    [1e7, 200],
    [1e8, 200],
]

params = fill(zeros(Int64, 3), length(dhw_scenarios) * length(n_seed_locations))
idx = 1
for idx_seed_loc in 1:length(n_seed_locations)
    for idx_dhw_scen in 1:length(dhw_scenarios)
        params[idx] = [
            n_seed_locations[idx_seed_loc][1],
            n_seed_locations[idx_seed_loc][2],
            dhw_scenarios[idx_dhw_scen]
        ]
        global idx += 1
    end
end


pd_frequency::Int64 = 5
scens = []
N_seed_weights = (N_seed_TA=0.15, N_seed_CA=0.5, N_seed_CNA=0.0, N_seed_SM=0.35, N_seed_LM=0.0)
for param in params
    ADRIA.fix_factor!(dom;
        N_seed_TA=param[1] * N_seed_weights.N_seed_TA,
        N_seed_CA=param[1] * N_seed_weights.N_seed_CA,
        N_seed_CNA=param[1] * N_seed_weights.N_seed_CNA,
        N_seed_SM=param[1] * N_seed_weights.N_seed_SM,
        N_seed_LM=param[1] * N_seed_weights.N_seed_LM,
        dhw_scenario=param[3],
        min_iv_locations=param[2])
    scen = ADRIA.sample_options(dom, pd_frequency)
    push!(scens, scen)
end
scens = vcat(scens...)

# Run scenarios
#rs = ADRIA.run_scenarios(dom, scens[1:2,:], RCP)
rs = ADRIA.run_scenarios(dom, scens, RCP)

# ----------------------------------------------------------
# Load scenarios
#rs = ADRIA.load_results("./Outputs/ReefMod Engine__RCPs_26__2025-04-23_18_46_19_610")
#
#rss = [
#    ADRIA.load_results("./Outputs/ReefMod Engine__RCPs_26__2025-04-16_21_39_48_820"),
#    ADRIA.load_results("./Outputs/ReefMod Engine__RCPs_45__2025-04-16_19_07_38_604"),
#    ADRIA.load_results("./Outputs/ReefMod Engine__RCPs_70__2025-04-16_22_09_59_770")
#]
rss = [
    ADRIA.load_results("./Outputs/ReefMod Engine__RCPs_26__2025-04-23_18_46_19_610"),
    ADRIA.load_results("./Outputs/ReefMod Engine__RCPs_45__2025-04-24_12_40_27_234"),
    ADRIA.load_results("./Outputs/ReefMod Engine__RCPs_70__2025-04-25_10_35_04_772")
]

options = DataFrame(;
    option_name=String[],
    pathway_diversity=Float64[],
    N_seed=Int64[],
    dhw_scenario=Int64[],
    n_locations=Int64[],
    rcp=String[]
)
rcps = ["26", "45", "70"]
for (idx_rcp, rcp) in enumerate(rcps)
    rs = rss[idx_rcp]
    for param in params
        println(param)
        condition =
            scens.N_seed_CA .== param[1] .&& scens.min_iv_locations .== param[2] .&&
            scens.dhw_scenario .== param[3]
        idx_scens = findall(condition)
        tmp_options = ADRIA.analysis.pathway_diversity(rs, scens, idx_scens)
        tmp_options.N_seed = fill(param[1], size(tmp_options, 1))
        tmp_options.n_locations = fill(param[2], size(tmp_options, 1))
        tmp_options.dhw_scenario = fill(param[3], size(tmp_options, 1))
        tmp_options.rcp = fill(rcp, size(tmp_options, 1))
        global options = vcat(options, tmp_options)
    end
end

CSV.write(joinpath(pd_config["plot_output_path"], "pathway_diversity.csv"), options)
options = CSV.read(joinpath(pd_config["plot_output_path"], "pathway_diversity.csv"), DataFrame)


# RCP plot
option_fix_seed = options[
    options.N_seed .== 1e8 .&& options.n_locations .== 200,
    [:option_name, :pathway_diversity, :dhw_scenario, :rcp]
]

option_fix_seed = combine(groupby(option_fix_seed, [:option_name, :rcp])) do subdf
    diversity_values = subdf.pathway_diversity
    (
        mean_pd=mean(diversity_values),
        min_pd=minimum(diversity_values),
        max_pd=maximum(diversity_values)
    )
end

# Create mapping columns
unique_options = unique(options.option_name)
option_map = Dict(opt => i for (i, opt) in enumerate(unique_options))
option_fix_seed.option_idx = [option_map[opt] for opt in option_fix_seed.option_name]

unique_rcp = unique(option_fix_seed.rcp)
rcp_map = Dict(rcp => i for (i, rcp) in enumerate(unique_rcp))
option_fix_seed.rcp_idx = [rcp_map[rcp] for rcp in option_fix_seed.rcp]

# Plot RCP figure
fig = Figure(; size=(800, 600))
ax = Axis(fig[1, 1];
    xlabel="Option Name",
    ylabel="Pathway Diversity",
    yticks=(4.5:0.05:5.0),
    xticks=(1:length(unique_options), string.(unique_options)),
    limits=(nothing, (4.5, 5.0))
)
palette = Makie.current_default_theme().palette.color[]
barplot!(
    ax,
    option_fix_seed.option_idx,
    option_fix_seed.mean_pd;
    dodge=option_fix_seed.rcp_idx,
    color=[palette[i] for i in option_fix_seed.rcp_idx]
)
crossbar!(
    ax,
    option_fix_seed.option_idx,
    option_fix_seed.mean_pd,
    option_fix_seed.min_pd, option_fix_seed.max_pd;
    dodge=option_fix_seed.rcp_idx,
    color=(:gray50, 0.8)
)
elements = [PolyElement(; polycolor=palette[i]) for i in 1:length(unique_rcp)]
Legend(fig[1, 2], elements, string.(unique_rcp), "RCP")
save(joinpath(pd_config["plot_output_path"], "pathway_diversity_rcp.png"), fig)

# seed plot
option_fix_rcp = options[
    options.rcp .== "45",
    [:option_name, :pathway_diversity, :dhw_scenario, :N_seed, :n_locations]
]
option_fix_rcp.seed =
    string.(option_fix_rcp.N_seed) .* "_" .* string.(option_fix_rcp.n_locations)
option_fix_rcp = option_fix_rcp[
    option_fix_rcp.seed .∈ [["100_15", "1000000_200", "100000000_200", "10000000000_1000"]],
    [:option_name,
        :pathway_diversity,
        :dhw_scenario,
        :seed]
]

option_fix_rcp = combine(groupby(option_fix_rcp, [:option_name, :seed])) do subdf
    diversity_values = subdf.pathway_diversity
    (
        mean_pd=mean(diversity_values),
        min_pd=minimum(diversity_values),
        max_pd=maximum(diversity_values)
    )
end

# Create mapping columns
option_fix_rcp.option_idx = [option_map[opt] for opt in option_fix_rcp.option_name]

unique_seed = unique(option_fix_rcp.seed)
seed_map = Dict(seed => i for (i, seed) in enumerate(unique_seed))
option_fix_rcp.seed_idx = [seed_map[seed] for seed in option_fix_rcp.seed]

# Plot seed figure
fig = Figure(; size=(800, 600))
ax = Axis(fig[1, 1];
    xlabel="Option Name",
    ylabel="Pathway Diversity",
    yticks=(4.5:0.05:5.0),
    xticks=(1:length(unique_options), string.(unique_options)),
    limits=(nothing, (4.5, 5.0))
)
palette = Makie.current_default_theme().palette.color[]
barplot!(
    ax,
    option_fix_rcp.option_idx,
    option_fix_rcp.mean_pd;
    dodge=option_fix_rcp.seed_idx,
    color=[palette[i] for i in option_fix_rcp.seed_idx]
)
crossbar!(
    ax,
    option_fix_rcp.option_idx,
    option_fix_rcp.mean_pd,
    option_fix_rcp.min_pd, option_fix_rcp.max_pd;
    dodge=option_fix_rcp.seed_idx,
    color=(:gray50, 0.8)
)
elements = [PolyElement(; polycolor=palette[i]) for i in 1:length(unique_seed)]
Legend(
    fig[1, 2],
    elements,
    ["100", "1 million", "100 million", "10 billion"],
    "Number of seeds"
)
save(joinpath(pd_config["plot_output_path"], "pathway_diversity_seed.png"), fig)

# Cover by scenario
rs = rss[2]

s_tac = ADRIA.metrics.scenario_total_cover(rs)
idx_fix_param = findall(
    scens.N_seed_CA .== 1e10 .&& scens.min_iv_locations .== 1000 .&&
    scens.dhw_scenario .== 7
)
s_tac_clean = s_tac[
    scenarios=idx_fix_param,
    timesteps=1:Int64(
        rs.inputs.seed_year_start[1] + rs.inputs.seed_years[1]
    )
]

clusters = zeros(Int64, length(idx_fix_param))
for (option, idx) in option_map
    idx_scens = findall(
        option_ts -> option_ts[Int64(rs.inputs.seed_year_start[1])] == Symbol(option),
        scens[idx_fix_param, :].option_ts
    )
    clusters[idx_scens] .= idx
end

fig_opts = Dict(:size => (800, 400))
axis_opts = Dict(
    :ylabel => "Total absolute cover [m²]",
    :xlabel => "Timesteps [years]"
)
opts = Dict{Symbol,Any}(:summarize => true)

tsc_fig = ADRIA.viz.clustered_scenarios(
    s_tac_clean, clusters; opts=opts, fig_opts=fig_opts, axis_opts=axis_opts
)
save(joinpath(pd_config["plot_output_path"], "scenarios_tac.png"), tsc_fig)
