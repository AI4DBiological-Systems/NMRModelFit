
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


project_folder = "/home/roy/MEGAsync/outputs/NMR/align/Serine-BMRB-700-20mM-mod"
r = 1

function loadregion!(Bs, project_folder::String)

    load_path = joinpath(project_folder, "region_$(r).bson")
    dic = BSON.load(load_path)
    minf = dic[:minf]
    minx = dic[:minx]
    ret = dic[:ret]
    w = dic[:w]

    for n = 1:length(Bs)

        if !isempty(dic[:κs_βs])

            tmp = convert(Vector{Vector{Float64}}, dic[:κs_βs][n])
            for i = 1:length(Bs[n].ss_params.κs_β)
                Bs[n].ss_params.κs_β[i][:] = tmp[i][:]
            end
            Bs[n].ss_params.d[:] = convert(Vector{Float64}, dic[:ds][n])

            if !isempty(Bs[n].β_singlets)
                Bs[n].β_singlets[:] = convert(Vector{Float64}, dic[:β_singletss][n])
                Bs[n].d_singlets[:] = convert(Vector{Float64}, dic[:d_singletss][n])
            end
        end
    end

    return minf, minx, rets, w
end

minf, minx, rets, w = loadregion!(Bs, project_folder)

y_cost = y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]

q2 = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs; w = w)

U_cost_rad = U_cost .* (2*π)
cost = sum( abs2.(q2.(U_cost_rad) - y_cost) ) # 1.6576135224225783 for serine-mod.




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
    display_plot_flag = true,
    canvas_size = (1000, 400))

# check by plotting to see if we get the same plot back, and cost.

function convertΔω0toΔcs(y::T, fs::T, SW::T)::T where T
    return y*SW/(2*π*fs)
end

d_cs = convertΔω0toΔcs.(Bs[1].ss_params.d, fs, SW)
