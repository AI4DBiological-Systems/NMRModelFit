function importβd!(Bs::Vector{CompoundType{T, SpinSysParamsType1{T}}},
    κs_β_an::Vector{Vector{Vector{T}}}, d_an::Vector{Vector{T}},
    β_singlets_an::Vector{Vector{T}}, d_singlets_an::Vector{Vector{T}}) where T

    # for n = 1:length(Bs)
    #
    #     if !isempty(κs_β_an)
    #
    #         for i = 1:length(Bs[n].ss_params.κs_β)
    #             Bs[n].ss_params.κs_β[i][:] = κs_β_an[n][i][:]
    #         end
    #         Bs[n].ss_params.d[:] = d_an[n]
    #
    #         if !isempty(Bs[n].β_singlets)
    #             Bs[n].β_singlets[:] = β_singlets_an[n]
    #             Bs[n].d_singlets[:] = d_singlets_an[n]
    #         end
    #     end
    # end
    for n = 1:length(Bs)
        Bs[n].ss_params.d[:] = d_an[n]

        for k = 1:length(Bs[n].ss_params.κs_β)
            Bs[n].ss_params.κs_β[k][:] = κs_β_an[n][k]
        end

        Bs[n].β_singlets[:] = β_singlets_an[n]
        Bs[n].d_singlets[:] = d_singlets_an[n]
    end

    return nothing
end

function importβd!(Bs::Vector{CompoundType{T, SpinSysParamsType2{T}}},
    κs_β_an::Vector{Vector{Vector{T}}},
    κs_d_an::Vector{Vector{Vector{T}}},
    d_an::Vector{Vector{Vector{T}}},
    β_singlets_an::Vector{Vector{T}},
    d_singlets_an::Vector{Vector{T}}) where T

    for n = 1:length(Bs)

        # loop over non-singlet systems.
        for i = 1:length(Bs[n].ss_params.d)
            Bs[n].ss_params.d[i][:] = d_an[n][i][:]
        end

        for i = 1:length(Bs[n].ss_params.κs_β)
            @assert length(κs_β_an[n][i]) == length(κs_d_an[n][i])

            Bs[n].ss_params.κs_β[i][:] = κs_β_an[n][i]
            Bs[n].ss_params.κs_d[i][:] = κs_d_an[n][i]
        end

        Bs[n].β_singlets[:] = β_singlets_an[n]
        Bs[n].d_singlets[:] = d_singlets_an[n]
    end

    return nothing
end

"""
Update Bs with the fit results for many regions, then deepcopy the values in B.
The resultant nested arrays makes easier export/import/storage.
"""
function exportβdfromregionsresults!(Bs::Vector{CompoundType{T, SpinSysParamsType1{T}}},
    minx_ar::Vector{Vector{T}},
    obj_func_ar) where T

    N_regions = length(minx_ar)
    κs_β_an_ar = Vector{Vector{Vector{Vector{Float64}}}}(undef, N_regions)
    d_an_ar = Vector{Vector{Vector{Float64}}}(undef, N_regions)
    β_singlets_an_ar = Vector{Vector{Vector{Float64}}}(undef, N_regions)
    d_singlets_an_ar = Vector{Vector{Vector{Float64}}}(undef, N_regions)

    for r = 1:N_regions
        obj_func_ar[r](minx_ar[r])

        κs_β_an_ar[r] = collect( deepcopy(Bs[n].ss_params.κs_β) for n = 1:length(Bs))
        d_an_ar[r] = collect( deepcopy(Bs[n].ss_params.d) for n = 1:length(Bs))
        β_singlets_an_ar[r] = collect( deepcopy(Bs[n].β_singlets) for n = 1:length(Bs))
        d_singlets_an_ar[r] = collect( deepcopy(Bs[n].d_singlets) for n = 1:length(Bs))
    end

    return κs_β_an_ar, d_an_ar, β_singlets_an_ar, d_singlets_an_ar
end

"""
Update Bs with the fit results for many regions, then deepcopy the values in B.
The resultant nested arrays makes easier export/import/storage.
"""
function exportβdfromregionsresults!(Bs::Vector{CompoundType{T, SpinSysParamsType2{T}}},
    minx::Vector{T},
    obj_func) where T

    obj_func(minx)

    κs_β_an = collect( deepcopy(Bs[n].ss_params.κs_β) for n = 1:length(Bs))
    d_an = collect( deepcopy(Bs[n].ss_params.d) for n = 1:length(Bs))
    κs_d_an = collect( deepcopy(Bs[n].ss_params.κs_d) for n = 1:length(Bs))
    β_singlets_an = collect( deepcopy(Bs[n].β_singlets) for n = 1:length(Bs))
    d_singlets_an = collect( deepcopy(Bs[n].d_singlets) for n = 1:length(Bs))

    return κs_β_an, d_an, κs_d_an, β_singlets_an, d_singlets_an
end

#### conversion of cs to type2's κs_d.

# κs_d[n][i][l]
# cs_an[n][i][l]
function getcs(As::Vector{SHType{T}}, Phys::Vector{PhysicalParamsType{T}}) where T <: Real

    N_compounds = length(As)
    cs_an = Vector{Vector{Vector{T}}}(undef, N_compounds)
    cs_singlets_an = Vector{Vector{T}}(undef, N_compounds)

    for n = 1:N_compounds
        N_sys = length(As[n].Δc_bar)
        cs_an[n] = Vector{Vector{T}}(undef, N_sys)

        for i = 1:N_sys

            L = length(As[n].Δc_bar[i])
            inds_set = collect( l for l = 1:L )
            if !isempty(Phys[n].ME[i])
                L = length(Phys[n].ME[i])
                inds_set = collect( Phys[n].ME[i][l][1] for l = 1:L )
            end
            cs_an[n][i] = collect( Phys[n].cs_sys[i][inds_set[l]] for l = 1:L )
        end


        cs_singlets_an[n] = copy(Phys[1].cs_singlets)
    end

    return cs_an, cs_singlets_an
end

function getΔcstype2(P_y, cost_inds_set::Vector{Vector{Int}},
    cs_an::Vector{Vector{Vector{T}}},
    cs_singlets_an::Vector{Vector{T}},
    d_an_ar::Vector{Vector{Vector{T}}},
    d_singlets_an_ar::Vector{Vector{Vector{T}}},
    fs::T, SW::T) where T <: Real

    N_compounds = length(cs_an)

    # allocate output.
    A = Vector{Vector{Vector{Vector{T}}}}(undef, N_compounds)
    for n = 1:N_compounds
        N_sys = length(cs_an[n])

        A[n] = Vector{Vector{Vector{T}}}(undef, N_sys)

        for i = 1:N_sys

            Q = length(cs_an[n][i])
            A[n][i] = Vector{Vector{T}}(undef, Q)

            for l = 1:Q

                A[n][i][l] = Vector{T}(undef, 0)
            end
        end
    end

    S = Vector{Vector{Vector{T}}}(undef, N_compounds)
    for n = 1:N_compounds
        N_sys = length(cs_singlets_an[n])

        S[n] = Vector{Vector{Vector{T}}}(undef, N_sys)

        for i = 1:N_sys
            S[n][i] = Vector{T}(undef, 0)
        end
    end

    # update output.
    for r = 1:length(cost_inds_set)

        Pr = P_y[cost_inds_set[r]]
        min_Pr = minimum(Pr)
        max_Pr = maximum(Pr)

        for n = 1:N_compounds
            N_sys = length(cs_an[n])

            #A[n] = Vector{Vector{Vector{T}}}(undef, N_sys)
            for i = 1:N_sys

                Q = length(cs_an[n][i])
                #A[n][i] = Vector{Vector{T}}(undef, Q)

                for l = 1:Q

                    #A[n][i][l] = Vector{T}(undef, 0)
                    if min_Pr < cs_an[n][i][l] <= max_Pr
                        #println("$(min_Pr), $(cs_an[n][i][l]), $(max_Pr)")

                        # this region is active for this n,i.
                        push!(A[n][i][l], convertΔω0toΔcs(d_an_ar[r][n][i], fs, SW))
                    end
                end
            end

            for i = 1:length(cs_singlets_an[n])
                if min_Pr < cs_singlets_an[n][i] <= max_Pr
                    push!(S[n][i], convertΔω0toΔcs(d_singlets_an_ar[r][n][i], fs, SW))
                end
            end

        end
    end

    # average.
    B = Vector{Vector{Vector{T}}}(undef, N_compounds)
    V = Vector{Vector{T}}(undef, N_compounds)
    for n = 1:N_compounds
        N_sys = length(cs_an[n])

        B[n] = Vector{Vector{T}}(undef, N_sys)
        for i = 1:N_sys

            Q = length(cs_an[n][i])
            B[n][i] = Vector{T}(undef, Q)

            for l = 1:Q

                B[n][i][l] = zero(T)
                if !isempty(A[n][i][l])
                    B[n][i][l] = sum( A[n][i][l] )./ length(A[n][i][l])
                end
            end
        end

        V[n] = zeros(T, length(cs_singlets_an[n]))
        for i = 1:length(cs_singlets_an[n])
            if !isempty(S[n][i])
                V[n][i] = sum( S[n][i] )./ length(S[n][i])
            end
        end
    end

    return B, V, A, S
end

# the κs_β and κs_d variables have an ordering that corresponds to As[n].Δc_bar.
# They are reshuffled versions of the outputs of getcs() and getΔcstype2()
function initializeΔcstype1(Bs::Vector{CompoundType{T,SpinSysParamsType1{T}}}) where T
    N_compounds = length(Bs)
    Δcs_offset = Vector{Vector{T}}(undef, N_compounds)
    Δcs_offset_singlets = Vector{Vector{T}}(undef, N_compounds)

    for n = 1:N_compounds

        Δcs_offset[n] = zeros(T, length(Bs[n].ss_params.d))
        Δcs_offset_singlets[n] = zeros(T, length(Bs[n].d_singlets))
    end

    return Δcs_offset, Δcs_offset_singlets
end
