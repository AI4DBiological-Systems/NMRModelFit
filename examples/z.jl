# package up the sequence:
# -fit type1, all regions.
# - extract Δcs from fit result, use it to create Δcs_offset2
# - use Δcs_offset2 for fitting type2, entire mixture.

# - later, extract the cs values from type2 fit.
# - later, allow getAs to use supplied cs values from extracted cs values.

import NMRDataSetup
import NMRHamiltonian
import NMRSignalSimulator
import NMRSpecifyRegions

include("../src/NMRModelFit.jl")
import .NMRModelFit

import NLopt
import PlotlyJS
using Plots; plotly()

using LinearAlgebra
using FFTW

import BSON
import JSON

import Statistics
import Random
Random.seed!(25)

include("helper.jl")

SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs_reduce.json"
#SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs_low_intensity_threshold.json"

surrogate_config_path = "/home/roy/Documents/repo/NMRData/input/surrogate_configs/select_compounds_SH_configs.json"
fit_config_path = "/home/roy/Documents/repo/NMRData/input/fit_configs/align_700MHz_type1_select_compounds.json"

# get mapping from molecule names to their spin system info json files.
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")


### Methionine. # TODO diagnose.
experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/misc/bmse000915_methionine"
project_name = "Methionine-BMRB-600-100mM-bmse000915"
molecule_names = ["L-Methionine";]

# ### Glucose.
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/NRC/misc/glucose/Sep-25-2018"
# project_name = "Glucose-NRC-2018"
# molecule_names = ["alpha-D-Glucose"; "beta-D-Glucose";]

# ### Serine.
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-700-20mM/L-Serine"
# project_name = "Serine-BMRB-700-20mM-bmse000885"
# molecule_names = ["L-Serine";]


project_base_folder = "/home/roy/MEGAsync/outputs/NMR/align"
#project_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate" # old.
project_folder = joinpath(project_base_folder, project_name)
isdir(project_folder) || mkpath(project_folder)

y, U_y, P_y, As, Bs, Bs2, fs, SW, ν_0ppm, cost_inds_set, Δsys_cs,
u_min, u_max, cost_inds, Phys = setupfitsession(SH_config_path,
    surrogate_config_path,
    fit_config_path,
    H_params_path,
    dict_compound_to_filename,
    experiment_full_path,
    project_name,
    project_folder,
    molecule_names;
    solvent_ppm_guess = 4.7,
    solvent_window_ppm = 0.1,
    offset_ppm = 0.3,
    Δcs_max_scalar_default = 0.2,
    unique_cs_atol = 1e-6,
    prune_combo_Δc_bar_flag = true,
    region_min_dist = 0.1)

# new normalization.
Z = maximum( maximum(abs.(y[cost_inds_set[r]])) for r = 1:length(cost_inds_set) )
y = y ./ Z

a_setp, b_setp, minxs,
    rets = NMRModelFit.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)
#
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)
U = LinRange(u_min, u_max, 50000)
P = hz2ppmfunc.(U)

#@assert 1==2

### first fit.
Δcs_offset, Δcs_offset_singlets = NMRModelFit.initializeΔcstype1(Bs)

obj_funcs, minfs, minxs, rets, ws = runfitregions(y,
    U_y,
    P_y,
    P, U,
    As,
    Bs,
    fs,
    SW,
    Δsys_cs, a_setp, b_setp, cost_inds_set,
    u_min, u_max, Δcs_offset, Δcs_offset_singlets;
    maxeval = 50,
    save_BSON_flag = true,
    save_plot_flag = true,
    display_flag = true)

### prepare cs.
Δcs_an_ravg, Δcs_an_singlets_ravg, Δcs_an_r, Δcs_singlets_an_r = preparetype2Δcsoffset(As, Phys,
    fs, SW, ν_0ppm,
    P_y,
    cost_inds_set,
    project_folder,
    "results_regions_type1.bson")


Δsys_cs_refine = Δsys_cs .* 0.05
#Δsys_cs_refine = Δsys_cs .* 0.1
Δcs_offset2 = Δcs_an_ravg
Δcs_offset_singlets2 = Δcs_an_singlets_ravg

### 2nd fit, for the entire mixture.
y_cost = y[cost_inds]
P_cost = P_y[cost_inds]
U_cost = U_y[cost_inds]

obj_func, minf, minx, ret, w, extracted_Δcs = runfitmixture(y_cost, P_cost, U_cost,
    As, Bs2, fs, SW, Δsys_cs_refine,
    Δcs_offset2, Δcs_offset_singlets2, a_setp, b_setp, P, U,
    project_name, project_folder, "mixture_fit_", "mixture_fit.bson";
    maxeval = 50,
    save_BSON_flag = true,
    save_plot_flag = true,
    display_plot_flag = true,
    N_viz = 50000,
    display_reduction_factor = 1,
    display_threshold_factor = 0.001/10,
    canvas_size = (1000, 400))

dummy = 1
