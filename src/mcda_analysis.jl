#=
Run GBR-wide simulation with ADRIA using the identified MCDA methods.

Identify which MCDA method produces outcomes above unguided scenarios, which are used as a baseline.
=#

include("src/common.jl")
using Bootstrap

RCP = "45"
dom = ADRIA.load_domain(pd_config["domain_path"], RCP; calib_params_fn=pd_config["coral_param_path"])

ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "FogCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "MCCriteriaWeights").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "Coral").fieldname)
ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "GrowthAcceleration").fieldname)

N_seed_weights = (N_seed_TA=0.15, N_seed_CA=0.5, N_seed_CNA=0.0, N_seed_SM=0.35, N_seed_LM=0.0)
N_seed_total = 1e7
ADRIA.fix_factor!(
    dom;
    # Environmental
    wave_scenario=1,
    dhw_scenario=11,
    # Intervention parameters
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

n_samples = 256
cf_scens = ADRIA.sample_cf(dom, n_samples)
ug_scens = ADRIA.sample_unguided(dom, n_samples)

ADRIA.fix_factor!(dom; guided=1.0)
cocoso_scens = ADRIA.sample_guided(dom, n_samples)

ADRIA.fix_factor!(dom; guided=2.0)
mairca_scens = ADRIA.sample_guided(dom, n_samples)

ADRIA.fix_factor!(dom; guided=3.0)
moora_scens = ADRIA.sample_guided(dom, n_samples)

ADRIA.fix_factor!(dom; guided=4.0)
piv_scens = ADRIA.sample_guided(dom, n_samples)

ADRIA.fix_factor!(dom; guided=5.0)
vikor_scens = ADRIA.sample_guided(dom, n_samples)

scens = vcat(cf_scens, ug_scens, cocoso_scens, mairca_scens, moora_scens, piv_scens, vikor_scens)

rs = ADRIA.run_scenarios(dom, scens, "45")

# rs = ADRIA.run_scenarios(dom, scens, "45")
rs = ADRIA.load_results("./Outputs/ReefMod__RCPs_45__2024-03-16_23_51_50_743")
stac = ADRIA.metrics.scenario_total_cover(rs)
ADRIA.viz.scenarios(rs, stac)

# Plotting trajectories
f = ADRIA.viz.scenarios(rs, stac[:, scens.guided .== -1.0])
f0 = ADRIA.viz.scenarios(rs, stac[:, scens.guided .== 0.0])
f1 = ADRIA.viz.scenarios(rs, stac[:, scens.guided .== 1.0])
f2 = ADRIA.viz.scenarios(rs, stac[:, scens.guided .== 2.0])
f3 = ADRIA.viz.scenarios(rs, stac[:, scens.guided .== 3.0])
f4 = ADRIA.viz.scenarios(rs, stac[:, scens.guided .== 4.0])
f5 = ADRIA.viz.scenarios(rs, stac[:, scens.guided .== 5.0])

# Bootstrapped assessment against unguided scenarios
n_timesteps = size(stac, 1)

# Dims: Timestep, bootstrap mean/median, lower bound, upper bound, method
guide_delta_mean = zeros(n_timesteps, 3, 5)
guide_delta_median = zeros(n_timesteps, 3, 5)

ug_scen_idx = findall(scens.guided .== 0.0)
for guided_method in 1:5
    guided_scen_idx = findall(scens.guided .== guided_method)

    b_sample = zeros(1000, 2)

    # Bootstrap each timestep
    for n in axes(stac, 1)
        for s in axes(b_sample, 1)
            shuf_set = shuffle(guided_scen_idx)
            ug_scen_idx = shuffle(ug_scen_idx)
            delta = stac[n, shuf_set] .- stac[n, ug_scen_idx]

            b_sample[s, 1] = mean(delta)
            b_sample[s, 2] = median(delta)
        end

        # bootstrap mean/median for each shuffled sample
        guide_delta_mean[n, :, guided_method] .= collect(Iterators.flatten(confint(bootstrap(mean, b_sample[:, 1], BalancedSampling(100)), PercentileConfInt(0.95))))
        guide_delta_median[n, :, guided_method] .= collect(Iterators.flatten(confint(bootstrap(median, b_sample[:, 2], BalancedSampling(100)), PercentileConfInt(0.95))))
    end

    mean_count = count(guide_delta_mean[:, 1, guided_method] .> 0.0)
    median_count = count(guide_delta_median[:, 1, guided_method] .> 0.0)

    # Report the number of times the mean/median delta was > 0
    @info """
        $(string(ADRIA.mcda_methods()[guided_method])) - N timesteps where delta > 0:
        mean: $(mean_count) ($(round((mean_count / n_timesteps) * 100.0, digits=2))%)
        median: $(median_count) ($(round((median_count / n_timesteps) * 100.0, digits=2))%)
    """
end

# N timesteps is 79
# ┌ Info:     CocosoMethod - N timesteps where delta > 0:
# │     mean: 24 (30.38%)
# └     median: 33 (41.77%)
# ┌ Info:     MaircaMethod - N timesteps where delta > 0:
# │     mean: 38 (48.1%)
# └     median: 39 (49.37%)
# ┌ Info:     MooraMethod - N timesteps where delta > 0:
# │     mean: 31 (39.24%)
# └     median: 40 (50.63%)
# ┌ Info:     PIVMethod - N timesteps where delta > 0:
# │     mean: 52 (65.82%)
# └     median: 49 (62.03%)
# ┌ Info:     VikorMethod - N timesteps where delta > 0:
# │     mean: 26 (32.91%)
# └     median: 38 (48.1%)
