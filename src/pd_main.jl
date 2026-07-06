#=
Pathway diversity main simulation.
Main parameters varied:
- RCP
- dhw_scenario
- min_iv_locations
- Number of seeds
=#

include("src/common.jl")
using GraphMakie, CairoMakie, SankeyMakie

rcps = ["26", "45", "70"]
seed_years = 20
dom = ADRIA.load_domain(
    pd_config["domain_path"], rcps[1];
    calib_params_fn=pd_config["coral_param_path"],
    # timeframe: seed_years + 2 (start seeding), 5 (extra years)
    timeframe=(2022, 2022 + seed_years + 2 + 5)
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
    [1e8, 200]
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
N_seed_weights = (
    N_seed_TA=0.15, N_seed_CA=0.5, N_seed_CNA=0.0, N_seed_SM=0.35, N_seed_LM=0.0
)
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
rs = ADRIA.run_scenarios(dom, scens, rcps)

# ----------------------------------------------------------
# Load scenarios
# path = ""
# rs = ADRIA.load_results(path)

# Pathway diversity options ordering and decision metadata
option_names = ADRIA.analysis.option_seed_preference().option_name
n_options = length(option_names)
seed_year_start = Int64(rs.inputs.seed_year_start[1])

# All (param, rcp) combinations, in the same order used to fill the tables below
param_rcp = [(param, rcp) for param in params for rcp in rcps]

# Per-option pathway diversity table, pre-allocated with all params/rcp/option rows
n_rows = length(param_rcp) * n_options
options = DataFrame(;
    option_name=Vector{String}(undef, n_rows),
    pathway_diversity=zeros(n_rows),
    N_seed=Vector{Int64}(undef, n_rows),
    dhw_scenario=Vector{Int64}(undef, n_rows),
    n_locations=Vector{Int64}(undef, n_rows),
    rcp=Vector{String}(undef, n_rows)
)
for (block, (param, rcp)) in enumerate(param_rcp)
    base = (block - 1) * n_options
    for (i, option) in enumerate(option_names)
        options.option_name[base + i] = string(option)
        options.N_seed[base + i] = Int64(param[1])
        options.dhw_scenario[base + i] = Int64(param[3])
        options.n_locations[base + i] = Int64(param[2])
        options.rcp[base + i] = rcp
    end
end

# Per-scenario probability table (all scenarios), filled inside the loop
scenario_probs = DataFrame(;
    option_ts=rs.inputs.option_ts,
    N_seed=round.(Int64,
        rs.inputs.N_seed_TA .+ rs.inputs.N_seed_CA .+ rs.inputs.N_seed_CNA .+
        rs.inputs.N_seed_SM .+ rs.inputs.N_seed_LM),
    dhw_scenario=Int64.(rs.inputs.dhw_scenario),
    n_locations=Int64.(rs.inputs.min_iv_locations),
    rcp=string.(Int64.(rs.inputs.RCP)),
    probability=zeros(size(rs.inputs, 1))
)

for (block, (param, rcp)) in enumerate(param_rcp)
    @info "param $(param), RCP $(rcp)"
    condition =
        rs.inputs.N_seed_CA .== param[1] * N_seed_weights.N_seed_CA .&&
        rs.inputs.N_seed_TA .== param[1] * N_seed_weights.N_seed_TA .&&
        rs.inputs.min_iv_locations .== param[2] .&&
        rs.inputs.dhw_scenario .== param[3] .&&
        rs.inputs.RCP .== parse(Float64, rcp)
    idx_scens = findall(condition)

    scenario_result = ADRIA.analysis.pathway_diversity(rs, idx_scens; scenario_probabilities=true)

    # Fill the global per-scenario probability table at the returned scenario indices
    scenario_probs.probability[scenario_result.scenario_idx] = scenario_result.probability

    # Aggregate per starting option and compute pathway diversity
    base = (block - 1) * n_options
    for (i, option) in enumerate(option_names)
        mask = [ts[seed_year_start] == option for ts in scenario_result.decoded_ts]
        probs = scenario_result.probability[mask]
        options.pathway_diversity[base + i] = sum(ADRIA.analysis._entropy.(probs))
    end
end

CSV.write(joinpath(pd_config["plot_output_path"], "pathway_diversity.csv"), options)
CSV.write(
    joinpath(pd_config["plot_output_path"], "scenario_probabilities.csv"), scenario_probs
)

options = CSV.read(
    joinpath(pd_config["plot_output_path"], "pathway_diversity.csv"), DataFrame
)
scenario_probs = CSV.read(
    joinpath(pd_config["plot_output_path"], "scenario_probabilities.csv"), DataFrame
)

min_pd = floor(minimum(options.pathway_diversity) * 0.98; digits=1)
max_pd = ceil(maximum(options.pathway_diversity) * 1.02; digits=1)

# RCP plot
option_fix_seed = options[
    options.N_seed .== 1e7 .&& options.n_locations .== 200 .&& options.rcp .!= 26,
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

unique_rcp = sort(unique(option_fix_seed.rcp))
rcp_map = Dict(rcp => i for (i, rcp) in enumerate(unique_rcp))
option_fix_seed.rcp_idx = [rcp_map[rcp] for rcp in option_fix_seed.rcp]

# Within each RCP, rank options by descending pathway diversity so the
# highest-value option appears first (leftmost)
option_fix_seed = combine(groupby(option_fix_seed, :rcp_idx)) do subdf
    subdf = sort(subdf, :mean_pd; rev=true)
    subdf.dodge_idx = 1:nrow(subdf)
    subdf
end
sort!(option_fix_seed, [:rcp_idx, :dodge_idx])

# Plot RCP figure
fig = Figure(; size=(800, 300))
ax = Axis(fig[1, 1];
    xlabel="RCP",
    ylabel="Pathway Diversity",
    yticks=(min_pd:round((max_pd - min_pd) / 5; digits=2):max_pd),
    xticks=(1:length(unique_rcp), string.(unique_rcp)),
    limits=(nothing, (min_pd, max_pd))
)
palette = Makie.current_default_theme().palette.color[]
barplot!(
    ax,
    option_fix_seed.rcp_idx,
    option_fix_seed.mean_pd;
    dodge=option_fix_seed.dodge_idx,
    color=[palette[i] for i in option_fix_seed.option_idx]
)

n_dodge = length(unique_options)
dodge_offsets = [(i - (n_dodge + 1) / 2) * (0.8 / n_dodge) for i in 1:n_dodge]
dodged_x =
    option_fix_seed.rcp_idx .+ [dodge_offsets[i] for i in option_fix_seed.dodge_idx]
errorbars!(
    ax,
    dodged_x,
    option_fix_seed.mean_pd,
    option_fix_seed.mean_pd .- option_fix_seed.min_pd,
    option_fix_seed.max_pd .- option_fix_seed.mean_pd;
    color=:black,
    whiskerwidth=6
)
elements = [PolyElement(; polycolor=palette[i]) for i in 1:length(unique_options)]
Legend(fig[1, 2], elements, string.(unique_options), "Option")
save(joinpath(pd_config["plot_output_path"], "pathway_diversity_rcp.png"), fig)

# seed plot
option_fix_rcp = options[
    options.rcp .== 45,
    [:option_name, :pathway_diversity, :dhw_scenario, :N_seed, :n_locations]
]
option_fix_rcp.seed =
    string.(option_fix_rcp.N_seed) .* "_" .* string.(option_fix_rcp.n_locations)
option_fix_rcp = option_fix_rcp[
    option_fix_rcp.seed .∈ [["100_15", "1000000_200", "10000000_200", "100000000_200"]],
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
fig = Figure(; size=(800, 300))
ax = Axis(fig[1, 1];
    xlabel="Option Name",
    ylabel="Pathway Diversity",
    yticks=(min_pd:round((max_pd - min_pd) / 5; digits=2):max_pd),
    xticks=(1:length(unique_options), string.(unique_options)),
    limits=(nothing, (min_pd, max_pd))
)
palette = Makie.current_default_theme().palette.color[]
barplot!(
    ax,
    option_fix_rcp.option_idx,
    option_fix_rcp.mean_pd;
    dodge=option_fix_rcp.seed_idx,
    color=[palette[i] for i in option_fix_rcp.seed_idx]
)

n_dodge = length(unique_seed)
dodge_offsets = [(i - (n_dodge + 1) / 2) * (0.8 / n_dodge) for i in 1:n_dodge]
dodged_x = option_fix_rcp.option_idx .+ [dodge_offsets[i] for i in option_fix_rcp.seed_idx]
errorbars!(
    ax,
    dodged_x,
    option_fix_rcp.mean_pd,
    option_fix_rcp.mean_pd .- option_fix_rcp.min_pd,
    option_fix_rcp.max_pd .- option_fix_rcp.mean_pd;
    color=:black,
    whiskerwidth=6
)
elements = [PolyElement(; polycolor=palette[i]) for i in 1:length(unique_seed)]
Legend(
    fig[1, 2],
    elements,
    ["100", "1 million", "10 million", "100 million"],
    "Number of seeds"
)
save(joinpath(pd_config["plot_output_path"], "pathway_diversity_seed.png"), fig)

# Cover by scenario # TODO: change plot
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

# ----------------------------------------------------------
# Sankey diagram of option pathways
# Visualises, for a given parameter set, how probability mass flows between options
# across the `number_changes` decision points. Layer k = decision point k, node height =
# marginal probability of an option at that step, ribbon width = joint probability of the
# A→B transition between consecutive steps.

number_changes = seed_years ÷ pd_frequency
max_time = size(rs.seed_log, :timesteps)
# Timesteps of the decision points within the seeding window
decision_steps = [seed_year_start + (k - 1) * pd_frequency for k in 1:number_changes]

_sci(n) = (e = floor(Int, log10(n)); c = n / 10^e; isinteger(c) ? "$(round(Int,c))e$e" : "$(c)e$e")

# Starting option of a scenario = the option active at the first decision point
starting_option(option_ts) = ADRIA.analysis.decode_option_ts(
    option_ts, seed_year_start, seed_years, pd_frequency, max_time
)[seed_year_start]

# Draw a Sankey of option pathways into an existing axis `ax`. `df` should already be
# filtered to the scenarios to include (e.g. a single starting option).
function draw_sankey!(ax, df, option_names, number_changes)
    n_opt = length(option_names)
    opt_idx = Dict(o => i for (i, o) in enumerate(option_names))
    palette = Makie.current_default_theme().palette.color[]
    node(k, o) = (k - 1) * n_opt + o  # unique node id per (decision point, option)

    # Joint A→B transition mass between consecutive decision points
    link_mass = [zeros(n_opt, n_opt) for _ in 1:(number_changes - 1)]
    for row in eachrow(df)
        ts = ADRIA.analysis.decode_option_ts(
            row.option_ts, seed_year_start, seed_years, pd_frequency, max_time
        )
        path = ts[decision_steps]
        for k in 1:(number_changes - 1)
            link_mass[k][opt_idx[path[k]], opt_idx[path[k + 1]]] += row.probability
        end
    end

    # Build (source, target, weight) connections
    connections = Tuple{Int,Int,Float64}[]
    for k in 1:(number_changes - 1), a in 1:n_opt, b in 1:n_opt
        m = link_mass[k][a, b]
        m <= 0 && continue
        push!(connections, (node(k, a), node(k + 1, b), m))
    end

    # Compact node ids so SankeyMakie has no isolated (unconnected) nodes
    used = sort(unique(vcat([c[1] for c in connections], [c[2] for c in connections])))
    remap = Dict(id => i for (i, id) in enumerate(used))
    connections = [(remap[s], remap[t], w) for (s, t, w) in connections]
    nodelabels = [string(option_names[(id - 1) % n_opt + 1]) for id in used]
    nodecolors = [palette[(id - 1) % n_opt + 1] for id in used]

    # Force every column to follow `option_names` order (top to bottom), independent of
    # link weights. SankeyMakie needs the *full* pairwise constraint set per layer, and
    # `b => a` (for a < b) places option a above option b.
    forceorder = Pair{Int,Int}[]
    for k in 1:number_changes, a in 1:n_opt, b in (a + 1):n_opt
        (haskey(remap, node(k, a)) && haskey(remap, node(k, b))) || continue
        push!(forceorder, remap[node(k, b)] => remap[node(k, a)])
    end

    sankey!(
        ax, connections;
        nodelabels=nodelabels, nodecolor=nodecolors,
        linkcolor=SankeyMakie.SourceColor(0.4), forceorder=forceorder
    )
    hidedecorations!(ax)
    hidespines!(ax)
    return ax
end

# One file per dhw scenario, at a representative seed budget and RCP. Each file stacks
# 5 isolated Sankeys (one per starting option) as rows so the pathways never overlap.
sankey_N_seed = 1e6
sankey_rcp = 45

for dhw in dhw_scenarios
    df_sankey = scenario_probs[
        (scenario_probs.N_seed .== sankey_N_seed) .&
        (scenario_probs.dhw_scenario .== dhw) .&
        (scenario_probs.rcp .== sankey_rcp),
        :
    ]
    fig = Figure(; size=(1100, 300 * n_options))
    Label(
        fig[0, 1],
        "Option pathways — N_seed $(_sci(sankey_N_seed)), RCP $(sankey_rcp), dhw $(dhw)";
        fontsize=18, font=:bold
    )
    for (i, option) in enumerate(option_names)
        df_opt = df_sankey[[starting_option(ts) == option for ts in df_sankey.option_ts], :]
        ax = Axis(fig[i, 1]; title="Start: $(option)", titlealign=:left)
        if isempty(df_opt)
            hidedecorations!(ax)
            hidespines!(ax)
            continue
        end
        draw_sankey!(ax, df_opt, option_names, number_changes)
    end
    save(
        joinpath(
            pd_config["plot_output_path"],
            "sankey_dhw$(dhw)_rcp$(sankey_rcp)_seed$(_sci(sankey_N_seed)).png"
        ),
        fig
    )
end

# ----------------------------------------------------------
# Boxplots of scenario probability distributions per starting option
# One row per (N_seed, RCP) configuration, all at the same DHW scenario. Each box
# aggregates the probabilities of every scenario that shares a starting option.

boxplot_dhw = 5  # change to inspect a different DHW scenario
boxplot_configs = [
    (N_seed=10_000_000, n_locations=200, rcp=45),
    (N_seed=10_000_000, n_locations=200, rcp=70),
    (N_seed=1_000_000, n_locations=200, rcp=45)
]

n_opt = length(option_names)
opt_idx = Dict(o => i for (i, o) in enumerate(option_names))
palette = Makie.current_default_theme().palette.color[]

fig = Figure(; size=(800, 230 * length(boxplot_configs)))
for (row, cfg) in enumerate(boxplot_configs)
    df = scenario_probs[
        (scenario_probs.N_seed .== cfg.N_seed) .&
        (scenario_probs.n_locations .== cfg.n_locations) .&
        (scenario_probs.dhw_scenario .== boxplot_dhw) .&
        (scenario_probs.rcp .== cfg.rcp),
        :
    ]
    cats = [opt_idx[starting_option(ts)] for ts in df.option_ts]
    colors = [palette[c] for c in cats]
    ax = Axis(
        fig[row, 1];
        title="N_seed $(_sci(cfg.N_seed)), $(cfg.n_locations) locs, RCP $(cfg.rcp), dhw $(boxplot_dhw)",
        yticks=(1:n_opt, string.(option_names)),
        xlabel=(row == length(boxplot_configs) ? "Switching probability" : "")
    )
    boxplot!(ax, cats, df.probability; orientation=:horizontal, color=colors)
end
save(joinpath(pd_config["plot_output_path"], "probability_boxplots.png"), fig)

# ----------------------------------------------------------
# Lock-in score across decision steps
# Diagnoses how likely a decision is to *keep* the current option (e.g. past
# :heat_stress → new :heat_stress) rather than switch, as the seeding window
# progresses. For each transition k (between consecutive decision points), within each
# starting-option distribution (which sums to 1) we sum the probability mass of pathways
# that keep the same option, then average over starting options. Median over dhw with
# min/max whiskers, one colored line per (N_seed, RCP) configuration.

# Configurable parameter sets. Each entry is one colored line; n_locations fixed at 200.
lockin_configs = [
    (N_seed=1_000_000,   n_locations=200, rcp=45),
    (N_seed=1_000_000,   n_locations=200, rcp=70),
    (N_seed=10_000_000,  n_locations=200, rcp=45),
    (N_seed=10_000_000,  n_locations=200, rcp=70),
    (N_seed=100_000_000, n_locations=200, rcp=45)
]

# Long table: one row per (config, dhw, decision_step)
lockin = DataFrame(;
    N_seed=Int64[], n_locations=Int64[], rcp=Int64[],
    dhw_scenario=Int64[], decision_step=Int64[], keep_prob=Float64[]
)

for cfg in lockin_configs
    for dhw in dhw_scenarios
        df = scenario_probs[
            (scenario_probs.N_seed .== cfg.N_seed) .&
            (scenario_probs.n_locations .== cfg.n_locations) .&
            (scenario_probs.dhw_scenario .== dhw) .&
            (scenario_probs.rcp .== cfg.rcp),
            :
        ]
        isempty(df) && continue

        # Decode each scenario's option path at the decision points
        paths = [
            ADRIA.analysis.decode_option_ts(
                ts, seed_year_start, seed_years, pd_frequency, max_time
            )[decision_steps]
            for ts in df.option_ts
        ]

        # Per transition k, average over starting options of the diagonal (keep) mass
        for k in 1:(number_changes - 1)
            keep_per_option = Float64[]
            for option in option_names
                mask = [p[1] == option for p in paths]
                any(mask) || continue
                push!(keep_per_option,
                    sum(df.probability[mask][[p[k] == p[k + 1] for p in paths[mask]]]))
            end
            isempty(keep_per_option) && continue
            push!(lockin, (
                cfg.N_seed, cfg.n_locations, cfg.rcp, dhw, k, mean(keep_per_option)
            ))
        end
    end
end

CSV.write(joinpath(pd_config["plot_output_path"], "lockin_scores.csv"), lockin)

# Aggregate over dhw: median with min/max whiskers
lockin_agg = combine(
    groupby(lockin, [:N_seed, :n_locations, :rcp, :decision_step])
) do subdf
    (
        median_keep=median(subdf.keep_prob),
        min_keep=minimum(subdf.keep_prob),
        max_keep=maximum(subdf.keep_prob)
    )
end

# Plot: x = decision step, y = keep probability, one colored line per config
fig = Figure(; size=(800, 400))
ax = Axis(fig[1, 1];
    xlabel="Decision step",
    ylabel="P(keep current option)",
    xticks=(1:(number_changes - 1), string.(1:(number_changes - 1)))
)
palette = Makie.current_default_theme().palette.color[]
n_cfg = length(lockin_configs)
# Small per-series x-offset so overlapping whiskers stay legible
dodge_offsets = [(i - (n_cfg + 1) / 2) * 0.04 for i in 1:n_cfg]

for (i, cfg) in enumerate(lockin_configs)
    sub = sort(
        lockin_agg[
            (lockin_agg.N_seed .== cfg.N_seed) .&
            (lockin_agg.n_locations .== cfg.n_locations) .&
            (lockin_agg.rcp .== cfg.rcp),
            :
        ],
        :decision_step
    )
    isempty(sub) && continue
    xs = sub.decision_step .+ dodge_offsets[i]
    scatterlines!(ax, xs, sub.median_keep; color=palette[i])
    errorbars!(
        ax, xs, sub.median_keep,
        sub.median_keep .- sub.min_keep, sub.max_keep .- sub.median_keep;
        color=palette[i], whiskerwidth=6
    )
end

elements = [LineElement(; color=palette[i]) for i in 1:n_cfg]
labels = ["$(_sci(cfg.N_seed)) · RCP$(cfg.rcp)" for cfg in lockin_configs]
Legend(fig[1, 2], elements, labels, "Parameter set")
save(joinpath(pd_config["plot_output_path"], "lockin_scores.png"), fig)
