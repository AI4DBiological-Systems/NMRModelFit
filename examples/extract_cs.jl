
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

include("../src/NMRModelFit.jl")
import .NMRModelFit

include("helper.jl")

# run prep_full.jl first. Refines a region given the coarse fit.




minf_ar, minx_ar, ret_ar, w_ar, κs_β_an_ar, d_an_ar,
    β_singlets_an_ar, d_singlets_an_ar,
    dic = loadregionsresultstype1(project_folder, "results_regions_type1.bson")


# N_d = sum( NMRModelFit.getNdvars(Bs[n]) for n = 1:length(Bs) )
#
# r = 1
# NMRModelFit.importβd!(Bs, κs_β_an_ar[r], d_an_ar[r],
# β_singlets_an_ar[r], d_singlets_an_ar[r])
#
# ## want non-warped cs (not ω_cs) in ppm.
# extracted_Δcs = ones(N_d) .* Inf
# NMRModelFit.extractmixtured!(extracted_Δcs, Bs, As, 1, fs, SW)
#
# extracted_Δcs2 = NMRModelFit.prepareΔcsoffset(Bs, Bs2, As, fs, SW)
#
# ΩS0 = NMRModelFit.getΩS(As)
# ΩS0_ppm = NMRModelFit.getPs(ΩS0, hz2ppmfunc)
#
# # y_cost = y[cost_inds_set[r]]
# # P_cost = P_y[cost_inds_set[r]]
# # U_cost = U_y[cost_inds_set[r]]

active_flags_ar = collect( NMRModelFit.findcoveredspinsystems(ΩS0_ppm, P_y[cost_inds_set[r]]) for r = 1:length(cost_inds_set) )

κs_β_an, d_an, β_singlets_an, d_singlets_an,
    w = NMRModelFit.combineregionsresultstype1(As, hz2ppmfunc, P_y, cost_inds_set,
    w_ar, κs_β_an_ar, d_an_ar, β_singlets_an_ar, d_singlets_an_ar)

NMRModelFit.importβd!(Bs, κs_β_an, d_an, β_singlets_an, d_singlets_an)


### extract Δcs, and prepare Δcs_offset2 for type2 fit over whole region.
N_d = sum( NMRModelFit.getNdvars(Bs[n]) for n = 1:length(Bs) )

## want non-warped cs (not ω_cs) in ppm.
extracted_Δcs = ones(N_d) .* Inf
NMRModelFit.extractmixtured!(extracted_Δcs, Bs, As, 1, fs, SW)

extracted_Δcs2 = NMRModelFit.prepareΔcsoffset(Bs, Bs2, As, fs, SW)

# TODO I am here. merge region for whole compound fit.

r = 1


y_cost = y[cost_inds]
P_cost = P_y[cost_inds]
U_cost = U_y[cost_inds]

import PyPlot
PyPlot.plot(P_y, real.(y), label = "data")
PyPlot.plot(P_cost, y_cost, label = "cost positions", "x")

PyPlot.legend()
