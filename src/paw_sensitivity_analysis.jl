#=
PAW sensitivity analysis considering four main parameters:
- dhw_scenario
- min_iv_locations
- N_seed
- Seed weights
=#

using Revise
using Infiltrator
using Statistics, CSV, DataFrames
using ADRIA

RME_path = "/home/vitor/Code/ADRIA.jl/DataPackages/GBR_MCB_GBR_2026-03-30_v080/"
calib_path = "/home/vitor/Code/ADRIA.jl/DataPackages/calibrated_params.nc"
RCP = "45" # RCP 26, 45, 70

dom = ADRIA.load_domain(RME_path, RCP; calib_params_fn = calib_path)

ms = ADRIA.model_spec(dom)

# Fix four classes of parameters
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "FogCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "MCCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "Coral").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "GrowthAcceleration").fieldname)

# Fix manually some Environmental, Intervention and Depth parameters
ADRIA.fix_factor!(dom;
    # EnvironmentalLayer parameters
    wave_scenario=1,
    # Intervention parameters
    guided=1,                  # Use CoCoSo method
    N_mc_settlers=0,           # No moving corals
    seeding_devices_per_m2=5,  # Check if makes sense.
    fogging=0.0,               # No fogging
    SRM=0.0,                   # No shadding
    a_adapt=0.0,               # No assisted adaptation
    a_adapt_ref=0.0,           # No assisted adaptation
    seed_years=30,
    seed_deployment_freq=1,    # Lower bound of sampling distribution. Seed every year
    plan_horizon=20.0,         # Upper bound of sampling distribution
    seed_year_start=2,         # Start as soon as possible
    seed_strategy=1,           # Periodic seed deployment
    # decision.DepthThresholds parameters
    depth_min=2.0,             # Lower bound of distribution
    depth_offset=25.0,         # Upper bound of distribution
    # SeedCriteriaWeights parameters
    seed_wave_stress=0.0
)

# Set bounds for parameters: min_iv_locations
ADRIA.set_factor_bounds!(dom;
    min_iv_locations=(10, 1000),
    seed_heat_stress=(0.0,1.0),
    seed_in_connectivity=(0.0,1.0),
    seed_out_connectivity=(0.0,1.0),
    seed_depth=(0.0,1.0),
    seed_coral_cover=(0.0,1.0),
    seed_cluster_diversity=(0.0,1.0),
    seed_geographic_separation=(0.0,1.0),
    #seed_coral_diversity=(0.0,1.0)
)

num_samples = 1024
N_seed_values = [1e6, 1e7, 5e7, 1e8]

scens = []
for N_seed_value in N_seed_values
    # Fix all N_seed parameters so that sum is N_seed_value
    ADRIA.fix_factor!(dom;
        N_seed_TA=N_seed_value  ÷ 5,
        N_seed_CA=N_seed_value  ÷ 5,
        N_seed_CNA=N_seed_value ÷ 5,
        N_seed_SM=N_seed_value  ÷ 5,
        N_seed_LM=N_seed_value  ÷ 5
    )
    scen = ADRIA.sample(dom, num_samples ÷ length(N_seed_values))
    push!(scens, scen)
end
scens = vcat(scens...)

rs = ADRIA.run_scenarios(dom, scens, RCP)

# Load scenario
path = "Output/"
rs = ADRIA.load_results(path)

# Compute scenario metrics
s_tac = ADRIA.metrics.scenario_total_cover(rs)
mean_s_tac = vec(mean(s_tac, dims=1))

s_rsv = ADRIA.metrics.scenario_rsv(rs)
mean_s_rsv = vec(mean(s_rsv, dims=1))

s_juves = ADRIA.metrics.scenario_relative_juveniles(rs)
mean_s_juves = vec(mean(s_juves, dims=1))

s_even = ADRIA.metrics.scenario_evenness(rs)
mean_s_even = vec(mean(s_even, dims=1))

metrics = [mean_s_tac, mean_s_rsv, mean_s_juves, mean_s_even]

using GeoMakie, GraphMakie, WGLMakie

# Some shared options for the example plots below
fig_opts = Dict(:size => (1600, 800))
# Factors of Interest
opts = Dict(
    :factors => [:dhw_scenario, :N_seed_TA, :min_iv_locations, :seed_heat_stress,
        :seed_heat_stress, :seed_in_connectivity, :seed_out_connectivity,
        :seed_depth, :seed_coral_cover, :seed_cluster_diversity,
        :seed_geographic_separation]
)

for (idx, metric) in enumerate(metrics)
    si = ADRIA.sensitivity.pawn(rs, metric)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.pawn!(g, si; opts)
    Colorbar(g[1, 2]; colormap=:viridis, colorrange=(0, 1), label="Relative Sensitivity", height=Relative(0.65))
    save("pawn_si_$(idx).png", f)
end
