# package up the sequence:
# -fit type1, all regions.
# - extract Δcs from fit result, use it to create Δcs_offset2
# - use Δcs_offset2 for fitting type2, entire mixture.

# - later, extract the cs values from type2 fit.
# - later, allow getAs to use supplied cs values from extracted cs values.

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






Phys = getphysicalparameters(molecule_names, base_path, dict_compound_to_filename;
    unique_cs_atol = 1e-6)
