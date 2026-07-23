using Revise, Infiltrator
using Statistics, CSV, DataFrames, YAXArrays, TOML
using ADRIA

const ROOT_PATH = dirname(@__DIR__)
_config_path = joinpath(ROOT_PATH, "config.toml")
_config = TOML.parsefile(_config_path)
pd_config = _config["PathwayDiversityAnalysis"]

"""
    fix_common_parameters!(dom)

Fix the model factors that are held constant, to the same value, across every analysis
script in `src/`. Call this right after `ADRIA.load_domain` and before any script-specific
`fix_factor!` calls.

Fixed here:
- Four full component groups: `FogCriteriaWeights`, `MCCriteriaWeights`, `Coral`,
  `GrowthAcceleration`
- Environmental: `wave_scenario`, `cyclone_mortality_scenario`
- Intervention: `guided` (CoCoSo), `plan_horizon`, and the alternative interventions that are
  switched off (`N_mc_settlers`, `fogging`, `SRM`, `mcb_duration`, `mcb_albedo`)
- Seeding: `seed_year_start`, `seed_deployment_freq`, `seed_strategy`, `a_adapt_ref`
- Depth thresholds: `depth_min`, `depth_offset`

Intentionally left to each script: factors that vary between scripts or in value, e.g.
`seed_years`, `dhw_scenario`, `min_iv_locations`, the `N_seed_*` budget, the per-option seed
weights, and — where they are varied (`pawn_sensitivity_analysis_counterfactual.jl`) —
`seeding_devices_per_m2` and `a_adapt`.
"""
function fix_common_parameters!(dom)
    ms = ADRIA.model_spec(dom)

    ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "FogCriteriaWeights").fieldname)
    ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "MCCriteriaWeights").fieldname)
    ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "Coral").fieldname)
    ADRIA.fix_factor!(dom, ADRIA.component_params(ms, "GrowthAcceleration").fieldname)

    ADRIA.fix_factor!(dom;
        # EnvironmentalLayer
        wave_scenario=1,
        cyclone_mortality_scenario=0,
        # Intervention
        guided=1,               # CoCoSo
        plan_horizon=5.0,
        # Alternative interventions (all switched off)
        N_mc_settlers=0,
        fogging=0.0,
        SRM=0.0,
        mcb_duration=0,
        mcb_albedo=0,
        # Seeding
        seed_year_start=2,
        seed_deployment_freq=1,
        seed_strategy=1,        # Periodic deployment
        a_adapt_ref=5,
        # decision.DepthThresholds
        depth_min=2.0,
        depth_offset=25.0
    )

    return dom
end
