#=
Script showing results to support pathway diversity simulation.
Plot selected locations, probability matrix and its subcomponents and distribution of selected locations. All for
the first time-step of the analysis.
=#

include("src/common.jl")
using WGLMakie, GeoMakie, GraphMakie

RCP = "45" # RCP 26, 45, 70
dom = ADRIA.load_domain(
    pd_config["domain_path"], RCP; calib_params_fn=pd_config["coral_param_path"]
)

# Basic decision_matrix parameters
mcda_method = ADRIA.decision.mcda_methods()[1]
min_locs = 200 # varies between 5 and 20

# cover
sum_cover = vec(sum(dom.init_coral_cover; dims=1).data)

# connectivitiy
k_area_locs = ADRIA.loc_k_area(dom)
area_weighted_conn = dom.conn.data .* k_area_locs
conn_cache = similar(area_weighted_conn)
in_conn, out_conn, network = ADRIA.connectivity_strength(
    area_weighted_conn, sum_cover, conn_cache; out_method=ADRIA.eigenvector_centrality
)

# Diversity in cluster
diversity_scores = ADRIA.decision.cluster_diversity(dom.loc_data.CB_CALIB_GROUPS)
separation_scores = ADRIA.decision.geographic_separation(dom.loc_data.mean_to_neighbor)

# dhw and wave
plan_horizon = 5
α = 0.99
decay = α .^ (1:(plan_horizon + 1)) .^ 2
dhw_scenario = 7
dhw_scens = dom.dhw_scens[scenarios=dhw_scenario,
    mcb_durations=ADRIA.At(dom.dhw_scens.mcb_durations[1]),
    albedo=ADRIA.At(dom.dhw_scens.albedo[1])]
dhw_projection = ADRIA.weighted_projection(dhw_scens, 1, plan_horizon, decay, 75)

#Diversity
n_groups = dom.coral_growth.n_groups
n_sizes = dom.coral_growth.n_sizes
n_locs = ADRIA.n_locations(dom)
init_cover_3d = ADRIA._reshape_init_cover(
    dom.init_coral_cover.data, (n_sizes, n_groups, n_locs)
)
loc_taxa_cover = ADRIA.metrics.relative_loc_taxa_cover(
    reshape(init_cover_3d, (1, n_groups, n_sizes, n_locs))
)
evenness = ADRIA.metrics.coral_evenness(loc_taxa_cover.data)
diversity = ADRIA.metrics.coral_diversity(evenness.data)[timesteps=1].data

# Valid locations
depth_criteria = ADRIA.identify_within_depth_bounds(dom.loc_data.depth_med, 2.0, 25.0)
valid_seed_locs = depth_criteria

decision_matrix = ADRIA.decision_matrix(
    dom.loc_ids[valid_seed_locs],
    collect(fieldnames(ADRIA.SeedCriteriaWeights));
    seed_depth=dom.loc_data.depth_med[valid_seed_locs],
    seed_cluster_diversity=diversity_scores[valid_seed_locs],
    seed_geographic_separation=separation_scores[valid_seed_locs],
    seed_coral_cover=sum_cover[valid_seed_locs],
    seed_heat_stress=dhw_projection[valid_seed_locs],
    #seed_wave_stress=wave_projection[valid_seed_locs],
    seed_in_connectivity=in_conn[valid_seed_locs],
    seed_out_connectivity=out_conn[valid_seed_locs],
    seed_coral_diversity=diversity[valid_seed_locs]
)

options = ADRIA.analysis.option_seed_preference(; include_weights=true)
number_options = size(options, 1)

# Options definition plot
fig = ADRIA.viz.mcda_options(options)
save(joinpath(pd_config["plot_output_path"], "pd_figures", "option_weight_matrix.png"), fig)

# Fil selected_locations for each option
options.selected_locations = [Vector{String}() for _ in 1:number_options]
for row in eachrow(options)
    row.selected_locations = ADRIA.select_locations(
        row.preference, decision_matrix, mcda_method, collect(1:length(dom.loc_ids)),
        min_locs
    )
end

plottable = GeoMakie.to_multipoly(dom.loc_data[:, :geometry])
labels = ["A", "B", "C", "D", "E"]
f = Figure(; size=(3 * 600, 2 * 900))
for (idx, row) in enumerate(eachrow(options))
    ga = GeoAxis(
        f[(idx - 1) ÷ 3 + 1, mod(idx - 1, 3) + 1];
        dest="+proj=latlong +datum=WGS84",
        title="$(labels[idx])) $(row.option_name)",
        titlesize=22,
        aspect=DataAspect(),
        limits=((141.8, 153.7), (-25.2, -9.8)),
        xgridwidth=0.5,
        ygridwidth=0.5
    )
    poly!(ga, plottable)
    points = ADRIA.analysis._filter_locations(dom.loc_data, row.selected_locations)
    points = ADRIA.centroids(points)
    scatter!(ga, points; color=:red, markersize=5)
    Label(
        f[(idx - 1) ÷ 3 + 1, mod(idx - 1, 3) + 1, Bottom()],
        "Longitude";
        padding=(0, 0, 0, -50),
        fontsize=16
    )
    Label(
        f[(idx - 1) ÷ 3 + 1, mod(idx - 1, 3) + 1, Left()],
        "Latitude";
        rotation=π / 2,
        padding=(0, 40, 0, 0),
        fontsize=16
    )
end
save(joinpath(pd_config["plot_output_path"], "pd_figures", "selected_locations.png"), f)

# Plot of switching probability matrix
number_options = size(options, 1)
prob_matrix = zeros(number_options, number_options)
for source in 1:number_options
    probs = ADRIA.analysis.switching_probability(
        options[source, :option_name], decision_matrix, dom.loc_data, mcda_method, min_locs
    )
    prob_matrix[source, :] = 100 * probs.probability
end

# Switching probability
fig = Figure()
ax = Axis(
    fig[1, 1];
    xlabel="To",
    ylabel="From",
    xticks=(1:number_options, string.(options.option_name)),
    yticks=(1:number_options, string.(options.option_name)),
    xticklabelrotation=2.0 / π,
    yreversed=true
)
limits = (10.0, 30.0)
heatmap!(ax, transpose(prob_matrix); colorrange=limits)
Colorbar(fig[1, 2]; limits=limits)
text!(ax,
    string.(round.(vec(prob_matrix); digits=1));
    position=[Point2f(x, y) for x in 1:number_options for y in 1:number_options],
    align=(:center, :center),
    color=:white,
    fontsize=14
)
save(joinpath(pd_config["plot_output_path"], "pd_figures", "probability_matrix.png"), fig)

# Switching probability subcomponents
ports = ADRIA.analysis._ports()

distance_port = zeros(number_options, number_options)
dispersion = zeros(number_options, number_options)
similarity_matrix = zeros(number_options, number_options)

for row in 1:number_options, col in 1:number_options
    row_locs = ADRIA.analysis._filter_locations(
        dom.loc_data, options[row, :selected_locations]
    )
    col_locs = ADRIA.analysis._filter_locations(
        dom.loc_data, options[col, :selected_locations]
    )
    common_ids = intersect(row_locs.UNIQUE_ID, col_locs.UNIQUE_ID)

    unique_row_locs = row_locs[row_locs.UNIQUE_ID .∉ [common_ids], :]
    unique_col_locs = col_locs[col_locs.UNIQUE_ID .∉ [common_ids], :]

    distance_port[row, col] = ADRIA.analysis.distance_port_score(
        unique_col_locs, unique_row_locs, ports
    )
    dispersion[row, col] = ADRIA.analysis.dispersion_score(col_locs, row_locs)
    similarity_matrix[row, col] = ADRIA.analysis.option_similarity(
        unique_col_locs, unique_row_locs
    )
end

function _heatmap_panel!(ax, mat, colorrange=(0.0, 1.0), digits=2)
    heatmap!(ax, transpose(mat); colorrange=colorrange)
    return text!(ax,
        string.(round.(vec(mat); digits=digits));
        position=[Point2f(x, y) for x in 1:number_options for y in 1:number_options],
        align=(:center, :center),
        color=:white,
        fontsize=12
    )
end

# 4-panel heatmap: switching probability, option similarity, distance to port, dispersion.
# prob_matrix is normalized to [0,1] for the shared colorbar; text still shows %.
option_labels = string.(options.option_name)
fig = Figure(; size=(2 * 520, 2 * 450))
ax1 = Axis(fig[1, 1]; title="Switching probability (%)", xlabel="To", ylabel="From",
    xticks=(1:number_options, option_labels), yticks=(1:number_options, option_labels),
    xticklabelrotation=2.0 / π, yreversed=true)
heatmap!(ax1, transpose(prob_matrix ./ 100); colorrange=(0.0, 1.0))
text!(ax1,
    string.(round.(vec(prob_matrix); digits=1));
    position=[Point2f(x, y) for x in 1:number_options for y in 1:number_options],
    align=(:center, :center),
    color=:white,
    fontsize=12
)

ax2 = Axis(fig[1, 2]; title="Option similarity", xlabel="To",
    xticks=(1:number_options, option_labels), yticks=(1:number_options, option_labels),
    xticklabelrotation=2.0 / π, yreversed=true)
_heatmap_panel!(ax2, similarity_matrix)

ax3 = Axis(fig[2, 1]; title="Distance to port", xlabel="To", ylabel="From",
    xticks=(1:number_options, option_labels), yticks=(1:number_options, option_labels),
    xticklabelrotation=2.0 / π, yreversed=true)
_heatmap_panel!(ax3, distance_port)

ax4 = Axis(fig[2, 2]; title="Dispersion", xlabel="To",
    xticks=(1:number_options, option_labels), yticks=(1:number_options, option_labels),
    xticklabelrotation=2.0 / π, yreversed=true)
_heatmap_panel!(ax4, dispersion)

Colorbar(fig[1:2, 3]; limits=(0.0, 1.0))
save(joinpath(pd_config["plot_output_path"], "pd_figures", "option_scores_matrix.png"), fig)

# 1-D plot: mean over rows for each destination option (col)
fig = Figure()
ax = Axis(
    fig[1, 1];
    xlabel="Option",
    ylabel="Mean score",
    xticks=(1:number_options, option_labels),
    yticks=(0.0:0.2:1.0),
    limits=(nothing, (0.0, 1.0))
)
scatterlines!(
    ax, 1:number_options, vec(mean(dispersion; dims=1)); color=:red, label="Dispersion"
)
scatterlines!(
    ax,
    1:number_options,
    vec(mean(distance_port; dims=1));
    color=:blue,
    label="Distance to port"
)
scatterlines!(
    ax,
    1:number_options,
    vec(mean(similarity_matrix; dims=1));
    color=:green,
    label="Option similarity"
)
scatterlines!(
    ax,
    1:number_options,
    vec(mean(prob_matrix; dims=1)) ./ 100;
    color=:orange,
    label="Switching prob."
)
axislegend(ax; position=:ct)
save(joinpath(pd_config["plot_output_path"], "pd_figures", "option_scores.png"), fig)

# Pathway diversity plot
options.pathway_diversity = zeros(size(options, 1))
n_scenario_option = size(rs.inputs, 1) ÷ size(options, 1)
#probs = ADRIA.analysis.pathway_diversity(rs, scens)

fig = Figure(; resolution=(1600, 400))
for idx in 1:number_options
    prob = probs[((idx - 1) * n_scenario_option + 1):(idx * n_scenario_option)]
    options[idx, :pathway_diversity] = sum(ADRIA.analysis._entropy.(prob))

    ax1 = Axis(fig[1, idx]; xlabel=string(options[idx, :option_name]))
    hist!(prob; bins=20)
end
save(
    joinpath(pd_config["plot_output_path"], "pd_figures", "probabilities_histogram.png"),
    fig
)

fig = Figure()
ax = Axis(
    fig[1, 1];
    xlabel="Option",
    ylabel="Pathway Diversity",
    xticks=(1:number_options, string.(options.option_name)),
    yticks=(6.5:0.1:7.0),
    limits=(nothing, (6.5, 7.0))
)
barplot!(1:number_options, options.pathway_diversity)
save(joinpath(pd_config["plot_output_path"], "pd_figures", "pathway_diversity.png"), fig)

fig_s_tac = ADRIA.viz.scenarios(
    rs, s_tac; fig_opts=fig_opts, axis_opts=Dict(:ylabel => "Scenario Total Cover")
)
save(joinpath(pd_config["plot_output_path"], "pd_figures", "total_cover.png"), fig_s_tac)
