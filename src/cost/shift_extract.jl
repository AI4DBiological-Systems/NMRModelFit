# the inverse of the methods in shift_fit.jl.

# Given Bs, update p.
function extractmixturedwarp!(p_buffer::Vector{T},
    Bs::Vector{CompoundType{T,SST}},
    As,
    p::Vector{T},
    st_ind::Int,
    fs::T,
    SW::T,
    Δsys_cs::Vector{Vector{T}},
    Δcs_offset::Vector{T},
    itp_a,
    itp_b)::Int where {T <: Real, SST}


    j = extractmixtured!(p_buffer, Bs, As, st_ind, fs, SW)
    j1 = Δcstoparameter!(p, Bs, p_buffer, st_ind, itp_a, itp_b, Δsys_cs, Δcs_offset)

    @assert j1 == j # debug.

    return j
end

function convertΔω0toΔcs(y::T, fs::T, SW::T)::T where T
    return y*SW/(2*π*fs)
end

##### coarse shift extraction. SpinSysParamsType1.

function extractmixtured!(p::Vector{T},
    Bs::Vector{CompoundType{T,SpinSysParamsType1{T}}},
    As,
    st_ind::Int,
    fs::T,
    SW::T)::Int where T <: Real

    j = st_ind - 1

    for n = 1:length(Bs)

        N_spins_sys = length(Bs[n].ss_params.d)

        for i = 1:N_spins_sys
            j += 1

            p[j] = convertΔω0toΔcs(Bs[n].ss_params.d[i], fs, SW)
        end

        for i = 1:length(Bs[n].d_singlets)
            j += 1

            p[j] = convertΔω0toΔcs(Bs[n].d_singlets[i], fs, SW)
        end
    end

    return j
end

# updates p with a transformed version of x.
function Δcstoparameter!(x::Vector{T},
    Bs::Vector{CompoundType{T, SpinSysParamsType1{T}}},
    p::Vector{T},
    st_ind::Int,
    itp_a,
    itp_b,
    Δsys_cs::Vector{Vector{T}},
    Δcs_offset::Vector{T}) where T <: Real

    j = st_ind - 1

    for n = 1:length(Bs)

        N_spins_sys = length(Bs[n].ss_params.d)
        @assert length(Δsys_cs[n]) == N_spins_sys + length(Bs[n].d_singlets)

        #first.
        if N_spins_sys > 0
            j += 1

            x[j] = (p[j]-Δcs_offset[j])/Δsys_cs[n][1]

            # itp.
            target = convertcompactdomain(x[j], -one(T), one(T), zero(T), one(T))
            a = itp_a(target)
            b = itp_b(target)
        end

        for i = 2:length(Bs[n].ss_params.d)
            j += 1

            p2 = (p[j]-Δcs_offset[j])/Δsys_cs[n][i]
            x[j] = MonotoneMaps.evalinversecompositelogisticprobit(p2, a, b, -one(T), one(T))
        end

        for i = 1:length(Bs[n].d_singlets)
            j += 1

            x[j] = (p[j]-Δcs_offset[j])/Δsys_cs[n][i+N_spins_sys]
        end
    end

    return j
end

function extractmixtured!(p::Vector{T},
    Bs::Vector{CompoundType{T,SpinSysParamsType2{T}}},
    As,
    st_ind::Int,
    fs::T,
    SW::T)::Int where T <: Real

    j = st_ind - 1

    for n = 1:length(Bs)

        N_spins_sys = length(Bs[n].ss_params.κs_d)

        for i = 1:N_spins_sys
            for l = 1:length(Bs[n].ss_params.κs_d[i])
                j += 1

                p[j] = convertΔω0toΔcs(Bs[n].ss_params.κs_d[i][l], fs, SW)
            end
        end

        for i = 1:length(Bs[n].d_singlets)
            j += 1

            p[j] = convertΔω0toΔcs(Bs[n].d_singlets[i], fs, SW)
        end
    end

    return j
end

#### conversion from extracted cs from a fit of type1 to Δcs_offset for a fit of type2.

# copies p1 's values to their corresponding entries in p2, for the cs-related structure in Bs.
function initializeΔcsoffset!(p2::Vector{T},
    Bs::Vector{CompoundType{T,SpinSysParamsType2{T}}},
    p1::Vector{T},
    As,
    st_ind1::Int,
    st_ind2::Int) where T <: Real

    j2 = st_ind2 - 1
    j1 = st_ind1 - 1

    for n = 1:length(As)

        N_spins_sys = length(Bs[n].ss_params.κs_d)

        for i = 1:N_spins_sys
            j1 += 1

            for l = 1:length(Bs[n].ss_params.κs_d[i])
                j2 += 1

                p2[j2] = p1[j1]
            end

        end

        for i = 1:length(Bs[n].d_singlets)
            j2 += 1
            j1 += 1

            p2[j2] = p1[j1]
        end
    end

    return j1, j2
end

function prepareΔcsoffset(Bs::Vector{CompoundType{T, SpinSysParamsType1{T}}},
    Bs2::Vector{CompoundType{T, SpinSysParamsType2{T}}},
    As::Vector{SHType{T}}, fs::T, SW::T) where T <: Real

    N_d = sum( getNdvars(Bs[n]) for n = 1:length(Bs) )

    extracted_Δcs = ones(N_d) .* Inf # debug: initialize to Inf.
    extractmixtured!(extracted_Δcs, Bs, As, 1, fs, SW)

    N_d2 = sum( getNdvars(Bs2[n]) for n = 1:length(Bs2) )
    Δcs_offset2 = ones(N_d2) .* Inf
    j1, j2 = initializeΔcsoffset!(Δcs_offset2, Bs2, extracted_Δcs, As, 1, 1)

    @assert j1 == N_d
    @assert j2 == N_d2

    return Δcs_offset2
end
