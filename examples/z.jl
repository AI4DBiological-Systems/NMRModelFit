# fit a single region at a time.
# loads the result of coarse fit, which is fit_regions.jl.

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


r = 1
T = Float64
N_starts = 500
local_optim_algorithm = NLopt.LN_BOBYQA
xtol_rel = 1e-9
maxeval = 200 # 2 # 50
maxtime = Inf
β_optim_algorithm = :GN_DIRECT_L
w_lb_default = 1e-1
w_ub_default = 100.0
β_max_iters = 500 # 2 # 500
β_xtol_rel = 1e-9
β_ftol_rel = 1e-9
β_maxtime = Inf

optim_algorithm = NLopt.GN_ESCH #NLopt.LN_BOBYQA
max_iters = 1000
#max_iters = 500
ftol_rel = 1e-9



y_cost = y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]
LS_inds = 1:length(U_cost)

U_rad_cost = U_cost .* (2*π)

# setup inner optim over β.
q0, updatedfunc, getshiftfunc, N_vars_set,
run_optim, obj_func_β, E_BLS, w_BLS, b_BLS, updateβfunc, updatewfunc,
q_β = NMRModelFit.setupcostnesteddwarpw(Bs, As, fs, SW, LS_inds, U_rad_cost,
    y_cost, Δsys_cs, a_setp, b_setp;
    β_optim_algorithm = β_optim_algorithm,
    w_lb_default = w_lb_default,
    w_ub_default = w_ub_default,
    β_max_iters = β_max_iters,
    β_xtol_rel = β_xtol_rel,
    β_ftol_rel = β_ftol_rel,
    β_maxtime = β_maxtime)
#
# set up outer optim over shifts.
N_β = sum( NMRModelFit.getNβ(Bs[n]) for n = 1:length(Bs) )
p_β = zeros(T, N_β) # persistant buffer.

f = pp->NMRModelFit.costnesteddw(U_rad_cost, y_cost, updatedfunc, updatewfunc,
pp, Bs, run_optim, E_BLS, w_BLS, b_BLS, p_β)

df = xx->FiniteDiff.finite_difference_gradient(f, xx)

opt = NLopt.Opt(optim_algorithm, N_β)

N_d = sum( NMRModelFit.getNd(Bs[n]) for n = 1:length(Bs) )
p0 = -ones(N_d) .* 0.018
#p0 = -zeros(N_d)

shift_lb = -ones(N_d) .* 0.02
shift_ub = ones(N_d) .* 0.02

minf, minx, ret, N_evals = NMRModelFit.runNLopt!(opt,
    p0,
    f,
    df,
    shift_lb,
    shift_ub;
    max_iters = max_iters,
    #max_iters = 1, # debug.
    xtol_rel = xtol_rel,
    ftol_rel = ftol_rel,
    maxtime = maxtime)

cost_eval0 = norm(reinterpret(T, E_BLS)*w_BLS - b_BLS)^2
cost2 = f(minx)
cost_eval = norm(reinterpret(T, E_BLS)*w_BLS - b_BLS)^2
println("abs(cost_eval-minf) = ", abs(cost_eval-minf))
println("abs(cost2-minf) = ", abs(cost2-minf))
println()


# p_test = ones(N_d)
# f(p_test)
#
# p_test = zeros(N_d)
# updatedfunc(p_test)
# d = Bs[1].ss_params.d
#
# p_test = ones(N_d)
# updatedfunc(p_test)
# d = Bs[1].ss_params.d


d = Bs[1].ss_params.d

c_bar = As[1].Δc_bar
κs_d = Bs[1].ss_params.κs_d

NMRModelFit.convertΔcstoΔω0.(p_test, fs, SW)

tmp = collect( dot(minx, As[1].Δc_bar[1][k]) for k = 1:length(As[1].Δc_bar[1]) )
NMRModelFit.convertΔcstoΔω0.(tmp, fs, SW)

#@assert 1==23

#### visualize.
N_viz = 50000

U = LinRange(u_min, u_max, N_viz)
P = hz2ppmfunc.(U)
U_rad = U .* (2*π)

q2 = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs; w = w_BLS )
q_U = q2.(U_rad)

plotregion(P, U, q_U, P_y, y, P_cost, y_cost, 0.0, 1,
    "", "$(cost_eval)", "";
    save_plot_flag = false,
    display_plot_flag = true,
    canvas_size = (1000,400))

@assert 1==2

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
