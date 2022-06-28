# combination of prep_save.jl and align_prep.jl.
# Loads from NMR experiment instead of intermediate BSON file.

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

# TODO phenylalanine and another serine 700 MHz with normal serine entry.

##### global constants.

#


SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs_reduce.json"
#SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs_low_intensity_threshold.json"

surrogate_config_path = "/home/roy/Documents/repo/NMRData/input/surrogate_configs/select_compounds_SH_configs.json"
fit_config_path = "/home/roy/Documents/repo/NMRData/input/fit_configs/align_700MHz_type1_select_compounds.json"

# get mapping from molecule names to their spin system info json files.
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")

# ###
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-700-20mM/L-Serine"
# project_name = "Serine-BMRB-700-20mM"
# molecule_names = ["L-Serine";]


###
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-700-20mM/L-Serine"
# project_name = "Serine-BMRB-700-20mM-mod"
# molecule_names = ["L-Serine - mod";]
# dummy_SSFID = NMRSignalSimulator.SpinSysParamsType1(0.0)

# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-700-20mM/L-Serine"
# project_name = "Serine-BMRB-700-20mM"
# molecule_names = ["L-Serine";]
# dummy_SSFID = NMRSignalSimulator.SpinSysParamsType1(0.0)

experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-700-20mM/L-Serine"
project_name = "Serine-BMRB-700-20mM-fine"
molecule_names = ["L-Serine";]
dummy_SSFID = NMRSignalSimulator.SpinSysParamsType2(0.0)

# TODO: disable warp and Δsys_cs_used multiplier, and use hard bounds.


# specify the NMR experiment folder
#experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/misc/bmse000297_ethanol/"
#experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/NRC/misc/glucose/Sep-25-2018"
#experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-500-0.5mM/L-Serine"
#experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/NRC/NRC_4_amino_acid_mixture_Jan_2022/1"
#experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-700-20mM/L-Isoleucine"
#experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-700-20mM/L-Glutamine"
#experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-500-0.5mM/L-Leucine"
#experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/misc/bmse000795_2_DSS"
#experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/misc/gissmo_DSS"
#experiment_full_path = "/home/roy/MEGAsync/outputs/NMR/experiments/misc/bmse000900_L-Phenylalanine"
#experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/leucine-600-bmse000920_1"
#experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-700-20mM/ATP"

#experiment_full_path = "/home/roy/MEGAsync/outputs/NMR/experiments/misc/bmse000900_L-Phenylalanine"
#experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/glucose-600-bmse000855_1"

# specify where the calibration results should be saved for this experiment.
#project_name = "ethanol"
#project_name = "NRC-glucose-2018"

#project_name = "Serine-BMRB-700-20mM-2000entry"
#project_name = "Serine-BMRB-700-20mM-mod"
#project_name = "Serine-BMRB-500-0.5mM-mod" #project_name = "Serine-BMRB-500-0.5mM"

 #project_name = "Serine-glucose-NRC-Jan2022"
#project_name = "NRC-Jan2022-serine-glucose-dss"

#project_name = "Isoleucine-BMRB-700-20mM"
#project_name = "Glutamine-BMRB-700-20mM"
#project_name = "Leucine-BMRB-500-0.5mM" # skipped.
#project_name = "Leucine-BMRB-600-100mM"
#project_name = "bmse000795_2_DSS" # skipped.
#project_name = "gissmo_DSS"
#project_name = "Phenylalanine-GISSMO-600-bmse000900"
#project_name = "Leucine-BMRB-600" #
#project_name = "ATP-BMRB-700-20mM"

#project_name = "Phenylalanine-GISSMO-600-bmse000900-poster"
#project_name = "glucose-600-bmse000855_1-poster"

#molecule_names = ["Ethanol";]
#molecule_names = ["alpha-D-Glucose"; "beta-D-Glucose";]
#molecule_names = ["L-Serine - mod"; "alpha-D-Glucose"; "beta-D-Glucose";]
#molecule_names = ["D-Glucose - 2000 study";] #molecule_names = ["D-(+)-Glucose";]
#molecule_names = ["D-(+)-Glucose";]


#molecule_names = ["L-Serine - 2000 study";]
#molecule_names = ["L-Serine - mod";]

#molecule_names = ["L-Isoleucine";]
#molecule_names = ["L-Glutamine";]
#molecule_names = ["L-Leucine";]
#molecule_names = ["DSS";]
#molecule_names = ["L-Phenylalanine"]

project_base_folder = "/home/roy/MEGAsync/outputs/NMR/align"
#project_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate" # old.
project_folder = joinpath(project_base_folder, project_name)
isdir(project_folder) || mkpath(project_folder)


##### user input for loading data.

## for loading the experiment.
solvent_ppm_guess = 4.7
solvent_window_ppm = 0.1

# start the first entry at the frequency corresponding to offset_ppm.
# This is a method to deal with wrap-around frequency.
offset_ppm = 0.3

## for surrogate construction.


#w = [1.0;] # relative concentration.

Δcs_max_scalar_default = 0.2 # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.

unique_cs_atol = 1e-6 # for assessing chemical equivalence.
prune_combo_Δc_bar_flag = true

## for determining cost function positions
region_min_dist = 0.1 # in ppm.

#####


##### load experiment, normalize data.
s_t, S, hz2ppmfunc, ppm2hzfunc, ν_0ppm, fs, SW, α_0ppm, β_0ppm, λ_0ppm, Ω_0ppm,
    results_0ppm, dic, α_solvent, β_solvent, λ_solvent, Ω_solvent,
    results_solvent = NMRDataSetup.loadspectrum(experiment_full_path;
    solvent_ppm = solvent_ppm_guess,
    solvent_window_ppm = solvent_window_ppm)
#

# start the first entry at the frequency corresponding to offset_Hz, which we set to 0.3 ppm.

offset_Hz = ν_0ppm - (ppm2hzfunc(offset_ppm)-ppm2hzfunc(0.0))
DFT_s = fft(s_t)
U_DFT, U_y, U_inds = NMRDataSetup.getwraparoundDFTfreqs(length(s_t), fs, offset_Hz)

S_U = DFT_s[U_inds]
P_y = hz2ppmfunc.(U_y)

val, ind = findmin( abs.(U_y .- ν_0ppm) )
Z = abs(S_U[ind])

y = S_U ./ Z
#####
#@assert 1==2

####### compute surrogate.
λ0 = λ_0ppm

# get a surrogate where K_{n,i} is encouraged to be no larger than `early_exit_part_size`.
println("Timing: getmageqinfomixture")
@time MEs = NMRHamiltonian.getmageqinfomixture(molecule_names,
    H_params_path,
    dict_compound_to_filename;
    unique_cs_atol = unique_cs_atol)

println("Timing: setupmixtureSH()")
@time mixture_params = NMRHamiltonian.setupmixtureSH(molecule_names,
    H_params_path, dict_compound_to_filename, fs, SW,
    ν_0ppm;
    MEs = MEs,
    config_path = SH_config_path,
    prune_combo_Δc_bar_flag = prune_combo_Δc_bar_flag)
As = mixture_params


ΩS_ppm = NMRModelFit.getΩSppm(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(NMRModelFit.combinevectors(ΩS_ppm))

u_offset = 0.5 # in ppm.
u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)



println("Timing: fitproxies!()")
Bs = NMRSignalSimulator.fitproxies(As, dummy_SSFID, λ0;
    names = molecule_names,
    config_path = surrogate_config_path,
    Δcs_max_scalar_default = Δcs_max_scalar_default,
    u_min = u_min,
    u_max = u_max)

####### end mixture proxy.


##### prepare fit positions.
### prepare positions.
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)
P_y = hz2ppmfunc.(U_y)

ΩS0 = NMRModelFit.getΩS(As)
ΩS0_ppm = NMRModelFit.getPs(ΩS0, hz2ppmfunc)

Δsys_cs, y_cost_all, U_cost_all, P_cost_all, exp_info, cost_inds, cost_inds_set,
    λ_lbs, λ_ubs, κs_β_orderings,
    κs_β_DOFs = NMRModelFit.prepareoptim(fit_config_path, molecule_names,
    hz2ppmfunc, U_y, y, As; region_min_dist = region_min_dist)
#
a_setp, b_setp, minxs,
    rets = NMRModelFit.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)

# new normalization.
Z = maximum( maximum(abs.(y[cost_inds_set[r]])) for r = 1:length(cost_inds_set) )
y = y ./ Z

# include("inner_kappa.jl") # debug inner optim.
#@assert 1==3

include("fit_regions.jl")
