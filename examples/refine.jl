# fit type1-regions, get Δcs_offset2 from fit results, then fit type2-single region or entire mixture.

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

save_plot_flag = true
display_plot_flag = true

# # no prior info.
# Δsys_cs_refine = Δsys_cs #.* 0.1
# Δcs_offset2 = zeros(length(extracted_Δcs2) )#-extracted_Δcs2

# uses prior optimization from type1.
Δsys_cs_refine = Δsys_cs .* 0.05
Δcs_offset2 = -extracted_Δcs2

N_d2 = sum( NMRModelFit.getNdvars(Bs2[n]) for n = 1:length(Bs2) )

shift_lb = -ones(N_d2)
shift_ub = ones(N_d2)

a_setp, b_setp, minxs,
    rets = NMRModelFit.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)

import NLopt

println("fitmodel")
@time obj_func, minf, minx, ret,
w = NMRModelFit.fitmodel(y_cost,
    U_cost,
    P_cost,
    As,
    Bs2,
    fs,
    SW,
    Δsys_cs_refine,
    a_setp, b_setp,
    shift_lb,
    shift_ub,
    Δcs_offset2;
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
    β_maxtime = Inf)

#
cost_obj_func = obj_func(minx)

# load these cs into SH simulation to get proper amplitude.
extracted_Δcs2 = ones(N_d2) .* Inf
NMRModelFit.extractmixtured!(extracted_Δcs2, Bs2, As, 1, fs, SW)

q2 = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs2; w = w)

U_cost_rad = U_cost .* (2*π)
cost = sum( abs2.(q2.(U_cost_rad) - y_cost) ) # 1.6576135224225783 for serine-mod.
println("abs(cost-cost_obj_func) = ", abs(cost-cost_obj_func))

#@assert 1==2

U = LinRange(u_min, u_max, 50000)
P = hz2ppmfunc.(U)
U_rad = U .* (2*π)

#ws[r][1] = 0.0
q_U = q2.(U_rad)

title_string = "$(project_name) type2 fit results, region $(r), real part, cost = $(cost)"
save_folder = project_folder
file_name = "type2_fit_real.html"

display_reduction_factor = 1
display_threshold_factor = 0.001/10
plotregion(P, U, q_U, P_y, y, P_cost, y_cost, display_threshold_factor, display_reduction_factor,
    save_folder, title_string, file_name;
    save_plot_flag = save_plot_flag,
    display_plot_flag = display_plot_flag,
    canvas_size = (1000, 400))
