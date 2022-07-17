
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

## manually modify.
#d_an[n][i][l]
#d_singlets_an[n][i]

# methionine.
d_singlets_an[1][1] = NMRModelFit.convertΔcstoΔω0(0.0004742469192480768, fs, SW)

d_an[1][1][5] = NMRModelFit.convertΔcstoΔω0(5.6805736901780626e-5, fs, SW) # region 1.

d_an[1][1][3] = NMRModelFit.convertΔcstoΔω0(0.025952319081648585, fs, SW) # 2.19 ppm region 3.
d_an[1][1][4] = NMRModelFit.convertΔcstoΔω0(0.025952319081648585, fs, SW) # 2.12 ppm region 3.

d_an[1][1][1] = NMRModelFit.convertΔcstoΔω0(0.0003895341441829677, fs, SW) # 3.8 ppm region 2.
d_an[1][1][2] = NMRModelFit.convertΔcstoΔω0(0.0003895341441829677, fs, SW) # 3.8 ppm region 2.


@time Phys11 = NMRHamiltonian.getphysicalparameters(["Ethanol"],
    H_params_path,
    dict_compound_to_filename;
    unique_cs_atol = 1e-6)

println("Timing: setupmixtureSH()")
@time As11 = NMRHamiltonian.setupmixtureSH(molecule_names,
    fs, SW, ν_0ppm,
    Phys;
    config_path = SH_config_path,
    prune_Δc_option = 5)

cs_an11, cs_singlets_an11 = NMRModelFit.getcs(As11, Phys11)

cs_an, cs_singlets_an = NMRModelFit.getcs(As, Phys)
Δcs_an_ravg, Δcs_an_singlets_ravg, Δcs_an_r, Δcs_singlets_an_r = NMRModelFit.getΔcstype2(P_y, cost_inds_set,
    cs_an, cs_singlets_an, d_an_ar, d_singlets_an_ar, fs, SW)

@assert 1==2

# region 2. 0.0003895341441829677
#  region 3. 0.025952319081648585  singlet 0.0004742469192480768
# region 1. 5.681598699376998e-5

# match which cs corresponds to which d_an.

# collect( NMRModelFit.convertΔω0toΔcs.(d_an[1][l],fs,SW) for l = 1:length(d_an) )
#
# NMRModelFit.convertΔω0toΔcs(0.2142,fs,SW)
# NMRModelFit.convertΔω0toΔcs(1.4688,fs,SW)
# NMRModelFit.convertΔω0toΔcs(97.859,fs,SW)

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

#@assert 1==2

minf_ar, minx_ar, ret_ar, w_ar, κs_β_an_ar, d_an_ar,
    β_singlets_an_ar, d_singlets_an_ar,
    dic = loadregionsresultstype1(project_folder, "results_regions_type1.bson")

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


### check cs.
N_d = sum( NMRModelFit.getNdvars(Bs[n]) for n = 1:length(Bs) )
extracted_Δcs = ones(N_d) .* Inf
NMRModelFit.extractmixtured!(extracted_Δcs, Bs, As, 1, fs, SW)

ΩS0 = NMRModelFit.getΩS(As)
ΩS0_ppm = NMRModelFit.getPs(ΩS0, hz2ppmfunc)
active_flags_ar = collect( NMRModelFit.findcoveredspinsystems(ΩS0_ppm, P_y[cost_inds_set[r]]) for r = 1:length(cost_inds_set) )
