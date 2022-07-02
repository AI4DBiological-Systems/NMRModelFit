
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

# run prep_full.jl first.

#project_folder = "/home/roy/MEGAsync/outputs/NMR/align/Serine-BMRB-700-20mM-mod"

r = 1


minf, minx, rets, w = loadregion!(Bs, project_folder)
κs_β = Bs[1].ss_params.κs_β
d = Bs[1].ss_params.d

N_d = sum( NMRModelFit.getNdvars(Bs[n]) for n = 1:length(Bs) )

y_cost = y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]

q2 = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs; w = w)

U_cost_rad = U_cost .* (2*π)
cost = sum( abs2.(q2.(U_cost_rad) - y_cost) ) # 1.6576135224225783 for serine-mod.




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
