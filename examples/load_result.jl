
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
#import PlotlyJS
#using Plots; plotly()

import BSON
import JSON

import Statistics
import Random
Random.seed!(25)

project_folder = "/home/roy/MEGAsync/outputs/NMR/align/Serine-BMRB-700-20mM-mod"
r = 1

load_path = joinpath(project_folder, "region_$(r).bson")
dic = BSON.load(load_path)
minf = dic[:minfs][r],
minx = dic[:minxs][r],
rets = dic[:rets][r],
w = dic[:w]

for n = 1:length(Bs)
    Bs[n].ss_params.κs_β = dic[:κs_βs][n]
    Bs[n].ss_params.d = dic[:ds][n]
    Bs[n].β_singlets = dic[:β_singletss][n]
    Bs[n].d_singlets = dic[:d_singletss][n]
end

# check by plotting to see if we get the same plot back, and cost.
