function importβd!(Bs::Vector{NMRModelFit.CompoundType{T, NMRModelFit.SpinSysParamsType1{T}}},
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

"""
Update Bs with the fit results for many regions, then deepcopy the values in B.
The resultant nested arrays makes easier export/import/storage.
"""
function exportβdfromregionsresults!(Bs::Vector{NMRModelFit.CompoundType{T, NMRModelFit.SpinSysParamsType1{T}}},
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
