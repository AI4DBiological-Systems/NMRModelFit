# fit a single region at a time.
# loads the result of coarse fit, which is fit_regions.jl.
# redo later.

include("../src/NMRModelFit.jl")
import .NMRModelFit

import NLopt
import PlotlyJS
using Plots; plotly()


include("helper.jl")


save_BSON_flag = true
save_plot_flag = true
# save_BSON_flag = false
# save_plot_flag = false

N_d = sum( NMRModelFit.getNd(Bs[n]) for n = 1:length(Bs) )
N_β = sum( NMRModelFit.getNβ(Bs[n]) for n = 1:length(Bs) )

shift_lb = -ones(N_d) .* 0.02
shift_ub = ones(N_d) .* 0.02

a_setp, b_setp, minxs,
    rets = NMRModelFit.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)
#

#LS_inds = 1:length(U_cost)
loop_range = 1:length(cost_inds_set)

println("Timing:")
@time obj_funcs, minfs, minxs, rets, ws = NMRModelFit.fitregions(y,
    U_y,
    P_y,
    As,
    Bs,
    fs,
    SW,
    Δsys_cs,
    a_setp, b_setp,
    shift_lb,
    shift_ub,
    cost_inds_set;
    loop_range = loop_range,
    N_starts = 100,
    local_optim_algorithm = NLopt.LN_BOBYQA,
    xtol_rel = 1e-9,
    maxeval = 50, # 2, # 50,
    maxtime = Inf,
    β_optim_algorithm = :GN_DIRECT_L,
    w_lb_default = 1e-1,
    w_ub_default = 100.0,
    β_max_iters = 500, # 2, # 500,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf);

dummy = 1

#minxs[3][1] = 0.0

#@assert 99==44

if save_BSON_flag
    save_path = joinpath(project_folder, "quantify_results.bson")
    BSON.bson(save_path, region_min_dist = region_min_dist,
    minfs = minfs,
    minxs = minxs,
    rets = rets,
    ws = ws)
end

#### visualize.
# minxs[1][1] = 0.000
# minxs[2][1] = 0.000
# display_reduction_factor = 100
# display_threshold_factor = 0.05/10
display_reduction_factor = 1
display_threshold_factor = 0.001/10
if "L-Isoleucine" in molecule_names
    display_reduction_factor = 1
    display_threshold_factor = 0.001/10
end
display_flag = true

# save data plot.
file_name = "data.html"
plots_save_path = joinpath(project_folder, file_name)



plotquantificationresults(As, Bs, ws, project_folder,
    P_y, y,
    region_min_dist,
    obj_funcs, minfs, minxs, rets,
    display_reduction_factor, display_threshold_factor,
    cost_inds_set, loop_range;
    canvas_size = (1000, 400),
    display_flag = display_flag,
    save_plot_flag = save_plot_flag,
    N_viz = 50000)

for r in loop_range
    println("region $(r):")
    println("objective: $(minfs[r]), return status: $(rets)")
    println("shift variable:")
    display(minxs)
    println()
end



# data plot.
file_name = "data.html"
plots_save_path = joinpath(project_folder, file_name)


canvas_size = (1000, 400)
plot_obj = Plots.plot( P_y,
    real.(y),
    title = "data spectrum",
    label = "data",
    seriestype = :line,
    ticks = :native,
    xlims = (P_y[1],P_y[end]),
    hover = P_y,
    linewidth = 4,
    xlabel = "ppm",
    ylabel = "real part of spectrum",
    xflip = true,
    size = canvas_size)

if save_plot_flag
    Plots.savefig(plot_obj, plots_save_path)
end
