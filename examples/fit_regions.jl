
include("../src/NMRModelFit.jl")
import .NMRModelFit

import NLopt
import PlotlyJS
using Plots; plotly()


include("helper.jl")


save_BSON_flag = true
save_plot_flag = true
display_flag = true

# save_BSON_flag = false
# save_plot_flag = false
# display_flag = true

Δsys_cs_used = deepcopy(Δsys_cs)

N_d = sum( NMRModelFit.getNdvars(Bs[n]) for n = 1:length(Bs) )
N_β = sum( NMRModelFit.getNβ(Bs[n]) for n = 1:length(Bs) )

Δcs_offset = zeros(N_d)

shift_lb = -ones(N_d)
shift_ub = ones(N_d)

a_setp, b_setp, minxs,
    rets = NMRModelFit.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)
#

#LS_inds = 1:length(U_cost)

println("Timing:")
@time obj_funcs, minfs, minxs, rets, ws = NMRModelFit.fitregions(y,
    U_y,
    P_y,
    As,
    Bs,
    fs,
    SW,
    Δsys_cs_used,
    a_setp, b_setp,
    shift_lb,
    shift_ub,
    cost_inds_set,
    Δcs_offset, Δcs_offset_singlets;
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




if save_BSON_flag
    save_name = "results_regions_type1.bson"
    saveregionsresults!(Bs, minfs, minxs, rets, ws, obj_funcs, project_folder, save_name)
end

#@assert 99==44

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


U = LinRange(u_min, u_max, N_viz)
P = hz2ppmfunc.(U)

plotquantificationresults(As, Bs, ws, project_folder,
    P_y, y, P, U,
    obj_funcs, minfs, minxs, rets,
    display_reduction_factor, display_threshold_factor,
    cost_inds_set;
    canvas_size = (1000, 400),
    display_flag = display_flag,
    save_plot_flag = save_plot_flag,
    N_viz = 50000,
    filename_prefix = "fit_type1_")

for r = 1:length(cost_inds_set)
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
