# alignquantificationregion
# costnesteddw

# assumes p is non-empty. x unused.
function applywarptoshifts!(x::Vector{T},
    Bs::Vector{CompoundType{T, SpinSysParamsType2{T}}},
    p::Vector{T},
    st_ind::Int,
    itp_a,
    itp_b) where T <: Real

    j = st_ind - 1

    for n = 1:length(Bs)

        ## for non-singlet spin systems.

        # set up warp function parameters, a, b.
        target = convertcompactdomain(p[j+1], -one(T), one(T), zero(T), one(T))
        a = itp_a(target)
        b = itp_b(target)

        # update.
        for i = 1:length(Bs[n].ss_params.κs_d)
            for l = 1:length(Bs[n].ss_params.κs_d[i])
                j += 1

                Bs[n].ss_params.κs_d[i][l] = MonotoneMaps.evalcompositelogisticprobit(p[j], a, b, -one(T), one(T))
            end
        end

        ## singlets.
        for i = 1:length(Bs[n].d_singlets)
            j += 1

            Bs[n].d_singlets[i] = p[j]
        end
    end

    return j
end

function updatemixtured!(Bs::Vector{CompoundType{T,SpinSysParamsType2{T}}},
    p::Vector{T}, # unused.
    As,
    st_ind::Int,
    fs::T,
    SW::T,
    Δsys_cs::Vector{Vector{T}})::Int where T <: Real

    #@assert length(warp_param_set) == length(Δ_shifts)

    j = st_ind - 1

    for n = 1:length(As)

        N_spins_sys = length(Bs[n].ss_params.d)
        @assert length(Δsys_cs[n]) == N_spins_sys + length(Bs[n].d_singlets)

        for i = 1:N_spins_sys

            for k = 1:length(Bs[n].ss_params.d)


                p2 = dot(Bs[n].ss_params.κs_β[i], As[n].Δc_bar[i][k])*Δsys_cs[n][i]
                Bs[n].ss_params.d[i][k] = convertΔcstoΔω0(p2, fs, SW)
            end
            j += length(Bs[n].ss_params.κs_β[i])
        end

        for i = 1:length(Bs[n].d_singlets)
            j += 1

            p2 = Bs[n].d_singlets[i]*Δsys_cs[n][i+N_spins_sys]
            Bs[n].d_singlets[i] = convertΔcstoΔω0(p2, fs, SW)
        end
    end

    return j
end
