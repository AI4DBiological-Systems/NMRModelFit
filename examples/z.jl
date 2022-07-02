# fit a single region at a time.
# loads the result of coarse fit, which is fit_regions.jl.

include("../src/NMRModelFit.jl")
import .NMRModelFit

import NLopt
import PlotlyJS
using Plots; plotly()


include("helper.jl")


####### tmp. create entry given new tag.
c2 = [3.97384; 3.93701; 3.83016]
c2m = c2 - [ -0.000586121; 0.00116044; 0.00137673]
c2p = c2 + [ -0.000586121; 0.00116044; 0.00137673]

c_mod = [3.9744; 3.9352; 3.83016]




function getphysicalparameters(target_names::Vector{String},
    base_path::String,
    dict_compound_to_filename;
    unique_cs_atol = 1e-6)

    N_compounds = length(target_names)
    H_IDss = Vector{Vector{Int}}(undef, N_compounds)
    Phys = Vector{NMRModelFit.PhysicalParamsType{Float64}}(undef, N_compounds)

    for n = 1:N_compounds

        # TODO add error-handling if name is not found in the dictionary, or filename does not exist.
        load_path = joinpath(base_path, dict_compound_to_filename[target_names[n]]["file name"])
        H_IDs, H_css, J_IDs, J_vals = NMRHamiltonian.loadcouplinginfojson(load_path)

        J_inds_sys, J_inds_sys_local, J_IDs_sys, J_vals_sys, H_inds_sys,
            cs_sys, H_inds_singlets, cs_singlets, H_inds, J_inds,
            g = NMRHamiltonian.setupcsJ(H_IDs, H_css, J_IDs, J_vals)

        ME = NMRHamiltonian.getmageqinfo(H_IDs, H_css, J_IDs, J_vals; unique_cs_atol = unique_cs_atol)

        Phys[n] = NMRModelFit.PhysicalParamsType(H_IDs, H_inds_sys, cs_sys,
            H_inds_singlets, cs_singlets, J_inds_sys, J_vals_sys, ME)
    end

    return Phys
end

Phys = getphysicalparameters(molecule_names, base_path, dict_compound_to_filename;
    unique_cs_atol = 1e-6)
