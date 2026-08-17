using CSV, DataFrames, CairoMakie

param = params[1]
rcp="45"
condition =
    rs.inputs.N_seed_CA .== param[1] * N_seed_weights.N_seed_CA .&&
    rs.inputs.N_seed_TA .== param[1] * N_seed_weights.N_seed_TA .&&
    rs.inputs.min_iv_locations .== param[2] .&&
    rs.inputs.dhw_scenario .== param[3] .&&
    rs.inputs.RCP .== parse(Float64, rcp)
idx_scens = findall(condition)
scenario_result = ADRIA.analysis.pathway_diversity(rs, idx_scens)

param = params[4]
rcp="45"
condition =
    rs.inputs.N_seed_CA .== param[1] * N_seed_weights.N_seed_CA .&&
    rs.inputs.N_seed_TA .== param[1] * N_seed_weights.N_seed_TA .&&
    rs.inputs.min_iv_locations .== param[2] .&&
    rs.inputs.dhw_scenario .== param[3] .&&
    rs.inputs.RCP .== parse(Float64, rcp)
idx_scens = findall(condition)
scenario_result = ADRIA.analysis.pathway_diversity(rs, idx_scens)

raw_rel_tac_low_seed = CSV.read("./Outputs/raw_rel_tac_seed_1e6.csv", DataFrame; header=false, types=Float64)
raw_rel_tac_high_seed = CSV.read("./Outputs/raw_rel_tac_seed_1e8.csv", DataFrame; header=false, types=Float64)

raw_fd_low_seed = CSV.read("./Outputs/raw_fd_seed_1e6.csv", DataFrame; header=false, types=Float64)
raw_fd_high_seed = CSV.read("./Outputs/raw_fd_seed_1e8.csv", DataFrame; header=false, types=Float64)

raw_rel_tac = (raw_rel_tac_low_seed, raw_rel_tac_low_seed_high_rcp, raw_rel_tac_high_seed)
raw_fd = (raw_fd_low_seed, raw_fd_low_seed_high_rcp, raw_fd_high_seed)

titles = ("seed=1e6, RCP=45", "seed=1e8, RCP=45")

fig_rel_tac = Figure(; size=(1200, 400))
for (i, data) in enumerate(raw_rel_tac)
    ax = Axis(fig_rel_tac[1, i]; title=titles[i], xlabel="raw", ylabel="count")
    hist!(ax, data.Column3; bins=40, color=(:steelblue, 0.7))
end
save("./Outputs/rel_tac_histogram.png", fig_rel_tac)

fig_fd = Figure(; size=(1200, 400))
for (i, data) in enumerate(raw_fd)
    ax = Axis(fig_fd[1, i]; title=titles[i], xlabel="raw", ylabel="count")
    hist!(ax, data.Column3; bins=40, color=(:darkorange, 0.7))
end
save("./Outputs/fd_histogram.png", fig_fd)


raw = -0.002:0.0001:0.002
scales = [0.001, 0.002, 0.003, 0.004, 0.005]
n = 0.2
f(x, scale) = 0.5 * (1 + tanh(x / scale))

fig = Figure()
ax = Axis(fig[1, 1]; xlabel="raw", ylabel="f(x)")
for scale in scales
    lines!(ax, raw, [f(x, scale) for x in raw]; label="scale = $scale")
end
axislegend(ax; position=:lt)
save("./Outputs/benefit_function.png", fig)




#=
# Log `raw` to a file depending on which call this is (identified by σ).
# Log to be placed inside two_sided_cvar. REMEMBER TO REMOVE THREADS
_raw_file = σ == _σ_rel_tac ? "./Outputs/raw_rel_tac.csv" : "./Outputs/raw_fd.csv"
open(_raw_file, "a") do io
    println(io, "$(cvar_lower),$(cvar_upper),$(raw)")
end
=#
