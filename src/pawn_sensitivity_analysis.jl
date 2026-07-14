#=
PAWN sensitivity analysis considering five main parameters:
- dhw_scenario
- min_iv_locations
- N_seed
- Seed weights (option)
- RCP
=#

include("src/common.jl")
using GeoMakie, GraphMakie, CairoMakie

seed_years = 20
rcps = ["26", "45", "70"]
dom = ADRIA.load_domain(
    pd_config["domain_path"], rcps[1];
    calib_params_fn=pd_config["coral_param_path"],
    # timeframe: seed_years + 2 (start seeding), 5 (extra years)
    timeframe=(2022, 2022 + seed_years + 2 + 5)
)
fix_common_parameters!(dom)

# Fix manually the remaining script-specific seeding parameters
ADRIA.fix_factor!(dom;
    seed_years=seed_years,
    seeding_devices_per_m2=5,
    a_adapt=5.0
)

N_seeds = [1e6, 1e7, 5e7, 1e8]
min_locations = [10, 50, 100, 500, 1000]
options = ADRIA.analysis.option_seed_preference(; include_weights=true)
dhw_scenarios = 1:11

n_scens = length(N_seeds) * nrow(options) * length(dhw_scenarios) * length(min_locations)
scens_base = repeat(ADRIA.sample(dom, 2)[1:1, :], n_scens)
scens_base[!, :option] = zeros(Int, n_scens)

row = 1
for N_seed in N_seeds, (opt_idx, option) in enumerate(eachrow(options)),
    dhw_scenario in dhw_scenarios, min_location in min_locations

    scens_base[row, :N_seed_TA] = N_seed ÷ 5
    scens_base[row, :N_seed_CA] = N_seed ÷ 5
    scens_base[row, :N_seed_CNA] = N_seed ÷ 5
    scens_base[row, :N_seed_SM] = N_seed ÷ 5
    scens_base[row, :N_seed_LM] = N_seed ÷ 5
    scens_base[row, :seed_heat_stress] = option[2]
    scens_base[row, :seed_in_connectivity] = option[3]
    scens_base[row, :seed_out_connectivity] = option[4]
    scens_base[row, :seed_depth] = option[5]
    scens_base[row, :seed_coral_cover] = option[6]
    scens_base[row, :seed_cluster_diversity] = option[7]
    scens_base[row, :seed_geographic_separation] = option[8]
    scens_base[row, :seed_coral_diversity] = option[9]
    scens_base[row, :dhw_scenario] = dhw_scenario
    scens_base[row, :min_iv_locations] = min_location
    scens_base[row, :option] = opt_idx
    row += 1
end

# Run one ResultSet per RCP, adding rcp_idx as a factor column, then combine
rs = ADRIA.run_scenarios(dom, scens_base, rcps)

# Load scenario
# path = "Output/"
# rs = ADRIA.load_results(path)

# Reefs that received seeding in every intervention scenario (across timesteps 2-32)
seed_per_reef_per_scen = dropdims(
    sum(
        rs.seed_log[timesteps=2:32, scenarios=1:nrow(rs.inputs)];
        dims=(:timesteps, :coral_id)
    );
    dims=(:timesteps, :coral_id)
)
always_seeded = vec(all(seed_per_reef_per_scen.data .> 0; dims=2))
never_seeded = vec(all(seed_per_reef_per_scen.data .== 0; dims=2))

selected_locations = .!(always_seeded .| never_seeded) # Remove reefs that always or never have seeding
selected_locations = findall(selected_locations)

s_tac = ADRIA.metrics.scenario_total_cover(rs; locations=selected_locations)
mean_s_tac = vec(mean(s_tac; dims=1))

s_rsv = ADRIA.metrics.scenario_rsv(rs; locations=selected_locations)
mean_s_rsv = vec(mean(s_rsv; dims=1))

s_even = ADRIA.metrics.scenario_evenness(rs; locations=selected_locations)
mean_s_even = vec(mean(s_even; dims=1))

# scenario_relative_juveniles ignores the locations kwarg, so pre-slice manually
_aj = ADRIA.metrics.absolute_juveniles(rs)
_k_area = ADRIA.loc_k_area(rs)[selected_locations]
s_juves = ADRIA.metrics.scenario_relative_juveniles(
    _aj[locations=selected_locations].data, _k_area
)
mean_s_juves = vec(mean(s_juves; dims=1))

metrics = [mean_s_tac, mean_s_rsv, mean_s_juves, mean_s_even]

# Some shared options for the example plots below
fig_opts = Dict(:size => (800, 400))
# Factors of Interest
opts = Dict(
    :factors => [:dhw_scenario, :RCP, :option, :min_iv_locations, :N_seed_TA],
    :by => :none,
    :ytick_labels => [
        "Climate Model",
        "RCP",
        "Option",
        "Number of locations",
        "Number of corals deployed"
    ]
)
axis_opts = Dict(
    :title => "",
    :titlesize => 22,
    :xticklabelsize => 20,
    :yticklabelsize => 20
)
titles = [
    "Total absolute cover",
    "Relative Shelter Volume",
    "Relative Juveniles",
    "Coral Evenness"
]

for (idx, metric) in enumerate(metrics)
    si = ADRIA.sensitivity.pawn(rs, metric)
    axis_opts[:title] = titles[idx]
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.pawn!(g, si; opts=opts, axis_opts=axis_opts)
    Colorbar(
        g[1, 2];
        colormap=:viridis,
        colorrange=(0, 1),
        label="Relative Sensitivity",
        labelsize=22,
        ticklabelsize=20,
        height=Relative(0.65)
    )
    title_clear = replace(lowercase(titles[idx]), " " => "_")
    save(joinpath(pd_config["plot_output_path"], "pawn_$(title_clear).png"), f)
end
