module NMRModelFit

#import JLD

using LinearAlgebra
import NLopt, BSON, Statistics, JSON, Interpolations

import MultistartOptimization

#import NMRData
import NMRDataSetup # https://github.com/AI4DBiological-Systems/NMRDataSetup.jl
import BoundedLeastSquares # https://github.com/nboyd/BoundedLeastSquares.jl

# dependencies to NMRSpecifyRegions. Need to add explicitly since these packages are not on the Julia public registry.
import NMRHamiltonian # https://github.com/AI4DBiological-Systems/NMRHamiltonian.jl
import NMRSignalSimulator # https://github.com/AI4DBiological-Systems/NMRSignalSimulator.jl

#
import NMRSpecifyRegions # https://github.com/AI4DBiological-Systems/NMRSpecifyRegions

import MonotoneMaps # https://github.com/RoyCCWang/MonotoneMaps.jl

include("../src/types.jl")
include("../src/utils.jl")

include("../src/front_end/coarse.jl")

include("../src/cost/prep.jl")
include("../src/cost/LSw.jl")
include("../src/cost/coarse_fit.jl")

end
