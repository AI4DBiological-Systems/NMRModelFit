
const SHType{T} = NMRHamiltonian.SHType{T} where T
const CompoundType{T,SST} = NMRSignalSimulator.CompoundType{T,SST} where {T,SST}
const SpinSysParamsType1{T} = NMRSignalSimulator.SpinSysParamsType1{T} where T
const SpinSysParamsType2{T} = NMRSignalSimulator.SpinSysParamsType2{T} where T


struct PhysicalParamsType{T}

    H_IDs::Vector{Int}

    H_inds_sys::Vector{Vector{Int}}
    cs_sys::Vector{Vector{Int}}

    H_inds_singlets::Vector{Int}
    cs_singlets::Vector{T}

    J_inds_sys::Vector{Vector{Tuple{Int,Int}}}
    J_vals_sys::Vector{Vector{T}}

    ME::Vector{Vector{Vector{Int}}}
end
