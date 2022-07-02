


# given p, update the corresponding attributes in Bs.
# assume shifts always between [-1,1]
# p_offset is the offset to p.
function updatemixturedwarp!(p_buffer::Vector{T},
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

    #p_buffer[:] = p_offset

    j1 = parametertoΔcs!(p_buffer, Bs, p, st_ind, itp_a, itp_b, Δsys_cs, Δcs_offset)
    j = updatemixtured!(Bs, p_buffer, As, st_ind, fs, SW)

    @assert j1 == j # debug.

    return j
end

function convertΔcstoΔω0(x::T, fs::T, SW::T)::T where T
    return x*2*π*fs/SW
end


###### coarse shift fit. SpinSysParamsType1.

# updates the Δcs container x with a transformed version of the parameter p.
function parametertoΔcs!(x::Vector{T},
    Bs::Vector{CompoundType{T, SpinSysParamsType1{T}}},
    p::Vector{T},
    st_ind::Int,
    itp_a,
    itp_b,
    Δsys_cs::Vector{Vector{T}},
    Δcs_offset::Vector{T}) where T <: Real

    j = st_ind - 1

    a = zero(T)
    b = zero(T)

    for n = 1:length(Bs)

        N_spins_sys = length(Bs[n].ss_params.d)
        @assert length(Δsys_cs[n]) == N_spins_sys + length(Bs[n].d_singlets)

        #first.
        if N_spins_sys > 0
            j += 1
            x[j] = p[j]*Δsys_cs[n][1] + Δcs_offset[j]

            # itp.
            target = convertcompactdomain(p[j], -one(T), one(T), zero(T), one(T))
            a = itp_a(target)
            b = itp_b(target)
        end

        for i = 2:N_spins_sys
            j += 1

            p2 = MonotoneMaps.evalcompositelogisticprobit(p[j], a, b, -one(T), one(T))
            x[j]  = p2*Δsys_cs[n][i] + Δcs_offset[j]
        end

        for i = 1:length(Bs[n].d_singlets)
            j += 1

            x[j] = p[j]*Δsys_cs[n][i+N_spins_sys] + Δcs_offset[j]
        end
    end

    return j
end

function updatemixtured!(Bs::Vector{CompoundType{T,SpinSysParamsType1{T}}},
    p::Vector{T},
    As,
    st_ind::Int,
    fs::T,
    SW::T)::Int where T <: Real

    #@assert length(warp_param_set) == length(Δ_shifts)

    j = st_ind - 1

    for n = 1:length(Bs)

        N_spins_sys = length(Bs[n].ss_params.d)

        for i = 1:N_spins_sys
            j += 1

            Bs[n].ss_params.d[i] = convertΔcstoΔω0(p[j], fs, SW)
        end

        for i = 1:length(Bs[n].d_singlets)
            j += 1

            Bs[n].d_singlets[i] = convertΔcstoΔω0(p[j], fs, SW)
        end
    end

    return j
end


################# fine shift fit. SpinSysParamsType2.

# assumes p is non-empty. x unused.
function parametertoΔcs!(x::Vector{T},
    Bs::Vector{CompoundType{T, SpinSysParamsType2{T}}},
    p::Vector{T},
    st_ind::Int,
    itp_a,
    itp_b,
    Δsys_cs::Vector{Vector{T}},
    Δcs_offset::Vector{T}) where T <: Real


    j = st_ind - 1

    for n = 1:length(Bs)

        N_spins_sys = length(Bs[n].ss_params.κs_d)
        @assert length(Δsys_cs[n]) == N_spins_sys + length(Bs[n].d_singlets)

        #first.
        if N_spins_sys > 0
            i = 1

            j += 1
            x[j] = p[j]*Δsys_cs[n][i] + Δcs_offset[j]

            # itp.
            target = convertcompactdomain(p[j], -one(T), one(T), zero(T), one(T))
            a = itp_a(target)
            b = itp_b(target)

            for l = 2:length(Bs[n].ss_params.κs_d[i])
                j += 1

                p2 = MonotoneMaps.evalcompositelogisticprobit(p[j], a, b, -one(T), one(T))
                x[j]  = p2*Δsys_cs[n][i] + Δcs_offset[j]
            end
        end

        for i = 2:N_spins_sys
            for l = 1:length(Bs[n].ss_params.κs_d[i])
                j += 1

                p2 = MonotoneMaps.evalcompositelogisticprobit(p[j], a, b, -one(T), one(T))
                x[j]  = p2*Δsys_cs[n][i] + Δcs_offset[j]
            end
        end

        ## singlets.
        for i = 1:length(Bs[n].d_singlets)
            j += 1

            x[j] = p[j]*Δsys_cs[n][i+N_spins_sys] + Δcs_offset[j]
        end
    end

    return j
end

# p is in units of chem shift.
function updatemixtured!(Bs::Vector{CompoundType{T,SpinSysParamsType2{T}}},
    p::Vector{T},
    As,
    st_ind::Int,
    fs::T,
    SW::T)::Int where T <: Real

    j = st_ind - 1

    for n = 1:length(As)

        N_spins_sys = length(Bs[n].ss_params.κs_d)

        for i = 1:N_spins_sys
            #for l = 1:length(Bs[n].ss_params.d[i])
            for l = 1:length(Bs[n].ss_params.κs_d[i])

                j += 1

                Bs[n].ss_params.κs_d[i][l] = p[j] # Δcs, not Δω0 = Δcs*2*π.
            end

            for k = 1:length(As[n].Δc_bar[i])

                # println("As[n].Δc_bar[$(i)][$(k)] = ", As[n].Δc_bar[i][k])
                p2 = dot(Bs[n].ss_params.κs_d[i], As[n].Δc_bar[i][k])#*Δsys_cs[n][i]
                Bs[n].ss_params.d[i][k] = convertΔcstoΔω0(p2, fs, SW)
            end
            #j += length(Bs[n].ss_params.κs_d[i])

            # println("Bs[n].ss_params.κs_d[$(i)] = ", Bs[n].ss_params.κs_d[i])
            # println()

        end

        for i = 1:length(Bs[n].d_singlets)
            j += 1

            Bs[n].d_singlets[i] = convertΔcstoΔω0(p[j], fs, SW)
        end
    end

    return j
end
