
import NMRDataSetup
import NMRHamiltonian
import NMRSignalSimulator
import NMRSpecifyRegions

include("../src/NMRModelFit.jl")
import .NMRModelFit
#import NMRModelFit

using LinearAlgebra
using FFTW

import PyPlot
import PlotlyJS
using Plots; plotly()

import BSON
import JSON

import Statistics
import Random
Random.seed!(25)

include("helper.jl")

# run prep_full.jl first. Refines a region given the coarse fit.

r = 1



minf, minx, rets, w = loadregion!(Bs, project_folder)
κs_β = Bs[1].ss_params.κs_β
d = Bs[1].ss_params.d

N_d = sum( NMRModelFit.getNd(Bs[n]) for n = 1:length(Bs) )

y_cost = y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]

# NMRModelFit.alignquantificationregion(y_cost,
#     U_cost,
#     P_cost,
#     As,
#     Bs,
#     fs,
#     SW,
#     Δsys_cs,
#     a_setp, b_setp,
#     shift_lb,
#     shift_ub,
#     Δcs_offset;
#     N_starts = N_starts,
#     local_optim_algorithm = NLopt.LN_BOBYQA,
#     xtol_rel = 1e-9,
#     maxeval = 50, # 2, # 50,
#     maxtime = Inf,
#     β_optim_algorithm = :GN_DIRECT_L,
#     w_lb_default = 1e-1,
#     w_ub_default = 100.0,
#     β_max_iters = 500, # 2, # 500,
#     β_xtol_rel = 1e-9,
#     β_ftol_rel = 1e-9,
#     β_maxtime = Inf)


q2 = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs; w = w)

U_cost_rad = U_cost .* (2*π)
cost = sum( abs2.(q2.(U_cost_rad) - y_cost) ) # 1.6576135224225783 for serine-mod.

a_setp, b_setp, minxs,
    rets = NMRModelFit.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)

Δcs_offset = zeros(N_d)

p_test = ones(N_d) .* Inf
x1 = ones(N_d) .* Inf
NMRModelFit.extractmixturedwarp!(x1, Bs, As, p_test, 1, fs, SW, Δsys_cs, Δcs_offset, a_setp, b_setp)


extracted_Δcs = ones(N_d) .* Inf
NMRModelFit.extractmixtured!(extracted_Δcs, Bs, As, 1, fs, SW)

# TODO I am here.

@assert 1==2

U = LinRange(u_min, u_max, 50000)
P = hz2ppmfunc.(U)
U_rad = U .* (2*π)

#ws[r][1] = 0.0
q_U = q2.(U_rad)

#file_name = "real_$(r).html"
title_string = "$(project_name) fit results, region $(r), real part, cost = $(cost)"
save_folder = ""
file_name = ""

display_reduction_factor = 1
display_threshold_factor = 0.001/10
plotregion(P, U, q_U, P_y, y, P_cost, y_cost, display_threshold_factor, display_reduction_factor,
    save_folder, title_string, file_name;
    save_plot_flag = false,
    display_plot_flag = false,
    canvas_size = (1000, 400))

# check by plotting to see if we get the same plot back, and cost.



d_cs = convertΔω0toΔcs.(Bs[1].ss_params.d, fs, SW)


T = Float64
x0 = 0.4
target = NMRModelFit.convertcompactdomain(x0, -one(T), one(T), zero(T), one(T))
a = a_setp(target)
b = b_setp(target)

warpfunc = tt->MonotoneMaps.evalcompositelogisticprobit(tt, a, b, -one(T), one(T))

t = LinRange(-1,1, 500)
y_t = warpfunc.(t)

PyPlot.plot(t, y_t)
