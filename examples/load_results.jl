
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


#project_folder = "/home/roy/MEGAsync/outputs/NMR/align/Serine-BMRB-700-20mM-mod"

### mixture results.
load_path = joinpath(project_folder, "mixture_fit.bson")
dic = BSON.load(load_path)

minf = convert(Float64, dic[:minf])
minx = convert(Vector{Float64}, dic[:minx])
ret = convert(Symbol, dic[:ret])
w = convert(Vector{Float64}, dic[:w])
κs_β_an = convert(Vector{Vector{Vector{Float64}}}, dic[:κs_β_an])
κs_d_an = convert(Vector{Vector{Vector{Float64}}}, dic[:κs_d_an])
d_an = convert(Vector{Vector{Vector{Float64}}}, dic[:d_an])
β_singlets_an = convert(Vector{Vector{Float64}}, dic[:β_singlets_an])
d_singlets_an = convert(Vector{Vector{Float64}}, dic[:d_singlets_an])


### region results.
minf_ar, minx_ar, ret_ar, w_ar, κs_β_an_ar, d_an_ar,
    β_singlets_an_ar, d_singlets_an_ar,
    dic = loadregionsresultstype1(project_folder, "results_regions_type1.bson")

cs_an, cs_singlets_an = NMRModelFit.getcs(As, Phys)
Δcs_an_ravg, Δcs_an_singlets_ravg, Δcs_an_r, Δcs_singlets_an_r = NMRModelFit.getΔcstype2(P_y, cost_inds_set,
    cs_an, cs_singlets_an, d_an_ar, d_singlets_an_ar, fs, SW)

#
Δcs_an_ar, Δcs_singlets_an_ar = NMRModelFit.getΔcstype1(d_an_ar, d_singlets_an_ar, fs, SW)

## manually modify mixture results.
#d_an[n][i][l]
#d_singlets_an[n][i]

# ##isoleucine.
# d_an[1][1][6] = NMRModelFit.convertΔcstoΔω0(0.002251401984308625, fs, SW) #.
# κs_β_an[1][1][:] = [-1.0457610662155323, 3.140156168608727, -1.0457610662155323, -1.0457610662155323, -1.0457610662155323, -2.4477704077352573]
# w[1] = 1.8075560511122795

# ## alpha-glucose. [n][i][k], last group is 5.2 ppm doublet.
# d_an[1][1][end] = NMRModelFit.convertΔcstoΔω0(-0.0056653618310591285, fs, SW)
#
# ## beta-glucose. [n][i][k] last group is 4.6 ppm doublet.
# d_an[2][1][end] = NMRModelFit.convertΔcstoΔω0(-0.006076603072532248, fs, SW)

# n = 2
# i = 1
# ordering, DOF = NMRHamiltonian.createorderingfromeqinds(Phys[n].ME[i], As[n].N_spins_sys[i])
# NMRModelFit.convertΔω0toΔcs.(d_an[2][1], fs, SW)


#@assert 1==2


# region 3.
#  region 2. 4.6
# region 1. 5.22


## update Bs2.
NMRModelFit.importβd!(Bs2, κs_β_an, κs_d_an, d_an, β_singlets_an, d_singlets_an)
q2 = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs2; w = w)


### plot.

y_cost = y[cost_inds]
P_cost = P_y[cost_inds]
U_cost = U_y[cost_inds]

U_cost_rad = U_cost .* (2*π)

U_rad = U .* (2*π)
q_U = q2.(U_rad)

cost2 = sum( abs2.(q2.(U_cost_rad) - y_cost) ) # 1.6576135224225783 for serine-mod.

display_reduction_factor = 1
display_threshold_factor = 0.001/10

file_name = ""
title_string = "mixture fit result, cost = $(cost2)"
save_folder = ""

plotregion(P, U, q_U, P_y, y, P_cost, y_cost,
    display_threshold_factor, display_reduction_factor,
    save_folder, title_string, file_name;
    save_plot_flag = false,
    display_plot_flag = true,
    canvas_size = (1000,400),
    spectrum_processing_func = real)

#
N_d2 = sum( NMRModelFit.getNdvars(Bs2[n]) for n = 1:length(Bs2) )
extracted_Δcs2 = ones(N_d2) .* Inf
NMRModelFit.extractmixtured!(extracted_Δcs2, Bs2, As, 1, fs, SW)


cs_an, cs_singlets_an = NMRModelFit.getcs(As, Phys)

Δcs_an_ravg, Δcs_an_singlets_ravg, Δcs_an_r, Δcs_singlets_an_r = preparetype2Δcsoffset(As, Phys,
    fs, SW, ν_0ppm,
    P_y,
    cost_inds_set,
    project_folder,
    "results_regions_type1.bson")

#@assert 21==2



r = 1
NMRModelFit.importβd!(Bs, κs_β_an_ar[r], d_an_ar[r], β_singlets_an_ar[r], d_singlets_an_ar[r])

q = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs; w = w_ar[r])

y_cost = y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]

U_cost_rad = U_cost .* (2*π)

U_rad = U .* (2*π)
q_U = q.(U_rad)

cost = sum( abs2.(q.(U_cost_rad) - y_cost) ) # 1.6576135224225783 for serine-mod.

display_reduction_factor = 1
display_threshold_factor = 0.001/10

file_name = ""
title_string = "region $(r) fit result, cost = $(cost)"
save_folder = ""

plotregion(P, U, q_U, P_y, y, P_cost, y_cost,
    display_threshold_factor, display_reduction_factor,
    save_folder, title_string, file_name;
    save_plot_flag = false,
    display_plot_flag = true,
    canvas_size = (1000,400),
    spectrum_processing_func = real)


# ### check cs.
# N_d = sum( NMRModelFit.getNdvars(Bs[n]) for n = 1:length(Bs) )
# extracted_Δcs = ones(N_d) .* Inf
# NMRModelFit.extractmixtured!(extracted_Δcs, Bs, As, 1, fs, SW)
#
# ΩS0 = NMRModelFit.getΩS(As)
# ΩS0_ppm = NMRModelFit.getPs(ΩS0, hz2ppmfunc)
# active_flags_ar = collect( NMRModelFit.findcoveredspinsystems(ΩS0_ppm, P_y[cost_inds_set[r]]) for r = 1:length(cost_inds_set) )
