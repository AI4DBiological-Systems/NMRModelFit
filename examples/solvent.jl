# fit solvent.

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

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])



##### global constants.

SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs_reduce.json"
#SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs_low_intensity_threshold.json"

surrogate_config_path = "/home/roy/Documents/repo/NMRData/input/surrogate_configs/select_compounds_SH_configs.json"
fit_config_path = "/home/roy/Documents/repo/NMRData/input/fit_configs/align_700MHz_type1_select_compounds.json"

# get mapping from molecule names to their spin system info json files.
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")

# ###

# ### Glucose.
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/NRC/misc/glucose/Sep-25-2018"
# project_name = "Glucose-NRC-2018"
# molecule_names = ["alpha-D-Glucose"; "beta-D-Glucose";]

### 4-amino acid.
experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/NRC/NRC_4_amino_acid_mixture_Jan_2022/1"
project_name = "NRC-Jan2022-serine-glucose-debug"
molecule_names = ["L-Serine"; "alpha-D-Glucose"; "beta-D-Glucose"]

# ### DMEM.
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/NRC/misc/dmem_medium/Oct-22-2012"
# project_name = "NRC-DMEM-solvent"
# molecule_names = ["L-Serine"; "alpha-D-Glucose"; "beta-D-Glucose"]




project_base_folder = "/home/roy/MEGAsync/outputs/NMR/align"
#project_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate" # old.
project_folder = joinpath(project_base_folder, project_name)
isdir(project_folder) || mkpath(project_folder)


##### user input for loading data.

## for loading the experiment.
solvent_ppm_guess = 4.7
solvent_window_ppm = 0.1

if project_name == "NRC-DMEM-solvent"
    solvent_ppm_guess = 4.9
end

# start the first entry at the frequency corresponding to offset_ppm.
# This is a method to deal with wrap-around frequency.
offset_ppm = 0.3

## for surrogate construction.


#w = [1.0;] # relative concentration.

??cs_max_scalar_default = 0.2 # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.

unique_cs_atol = 1e-6 # for assessing chemical equivalence.
prune_??c_option = 5

## for determining cost function positions
region_min_dist = 0.1 # in ppm.

#####


##### load experiment, normalize data.
s_t, S, hz2ppmfunc, ppm2hzfunc, ??_0ppm, fs, SW, ??_0ppm, ??_0ppm, ??_0ppm, ??_0ppm,
    results_0ppm, dic, ??_solvent, ??_solvent, ??_solvent, ??_solvent,
    results_solvent = NMRDataSetup.loadspectrum(experiment_full_path;
    solvent_ppm = solvent_ppm_guess,
    solvent_window_ppm = solvent_window_ppm)
#

# start the first entry at the frequency corresponding to offset_Hz, which we set to 0.3 ppm.

offset_Hz = ??_0ppm - (ppm2hzfunc(offset_ppm)-ppm2hzfunc(0.0))
DFT_s = fft(s_t)
U_DFT, U_y, U_inds = NMRDataSetup.getwraparoundDFTfreqs(length(s_t), fs, offset_Hz)

S_U = DFT_s[U_inds]
P_y = hz2ppmfunc.(U_y)

val, ind = findmin( abs.(U_y .- ??_0ppm) )
Z = abs(S_U[ind])

y = S_U ./ Z
#####

### solvent.
solvent_window_ppm = 0.3

??_solvent = ??_solvent ./ (2*??)
solvent_ppm = hz2ppmfunc(??_solvent)
??_lb = ppm2hzfunc(solvent_ppm - solvent_window_ppm)
??_ub = ppm2hzfunc(solvent_ppm + solvent_window_ppm)

??s_lb = [??_lb*2*??;]
??s_ub = [??_ub*2*??;]
??s_lb = [0.5*??_0ppm;]
??s_ub = [10*??_solvent;]

keep_flags = (??_lb .< U_y .< ??_ub)
U_cost_solvent = U_y[keep_flags]
y_cost_solvent = y[keep_flags]

??s_initial = [??_solvent;]
??s_initial = [??_solvent;]

U_cost_solvent_rad = U_cost_solvent .* (2*??)

minf, minx, ret, N_evals, f_LS, A_LS,
z_LS = NMRModelFit.fitsinglets(y_cost_solvent, U_cost_solvent_rad,
    ??s_lb, ??s_ub, ??s_lb, ??s_ub,
    ??s_initial, ??s_initial;
    max_iters = 2000,
    xtol_rel = 1e-9,
    ftol_rel = 0.0,#1e-9,
    maxtime = Inf)

#
f_LS(minx)
L = div(length(minx),2)
??s = minx[1:L]
??s = minx[L+1:end]

q_solvent = uu->NMRModelFit.evalclsinglets(uu, ??s, ??s, z_LS)
q_solvent_U = q_solvent.(U_cost_solvent_rad)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(hz2ppmfunc.(U_y), real.(y), label = "data")
PyPlot.plot(hz2ppmfunc.(U_cost_solvent), real.(q_solvent_U), label = "solvent")
PyPlot.plot(hz2ppmfunc.(U_cost_solvent), real.(q_solvent_U), "x")

PyPlot.gca().invert_xaxis()

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("solvent fit")


y2 = y - q_solvent.(U_y .* (2*??))
PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(hz2ppmfunc.(U_y), real.(y), label = "data")
PyPlot.plot(hz2ppmfunc.(U_y), real.(y2), label = "processed data")

PyPlot.gca().invert_xaxis()

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("data subtract solvent singlet")

# try DFT-FID singlet to see if the profile changes.

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(real.(s_t), label = "data")

PyPlot.legend()
PyPlot.xlabel("samples")
PyPlot.ylabel("real")
PyPlot.title("time domain data s_t")

# do Voigt fit for solvent.
# /home/roy/MEGAsync/work/literature/NMR/lineshape
# do DMEM z.jl

#### FID singlets.
function evalFIDsinglets(t_in::T, z::Vector{Complex{T}}, ??s, ??s, t0) where T <: Real

    t = t_in - t0

    out = zero(Complex{T})
    for i = 1:length(z)
        out += z[i]*cis(??s[i]*t)*exp(-??s[i]*t)
    end

    return out
end

t0 = 0.0
f = tt->evalFIDsinglets(tt, z_LS, ??s_star, ??s_star, t0)

@assert 1==4

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(real.(s_t), label = "data")

PyPlot.legend()
PyPlot.xlabel("samples")
PyPlot.ylabel("real")
PyPlot.title("time domain data s_t")
