## combine setup from multiple regions to one multi-interval region.
## combine results from multiple regions to one result.


# assumes st_a < fin_a, st_b < fin_b.
function isoverlapped(st_a::T, fin_a::T, st_b::T, fin_b::T) where T

    if st_a < st_b && st_b < fin_a < fin_b
        return true

    elseif st_b < st_a && st_a < fin_b < fin_a
        return true

    elseif st_b < st_a < fin_b
        return true

    elseif st_a < st_b < fin_a
        return true

    end

    return false
end

# Ωs0_ppm[n][i] contain the resonance freuqnecies for the n-th compound, i-th spin system.
# covered by P_cost.
function findcoveredspinsystems(ΩS0_ppm::Vector{Vector{Vector{T}}}, P_cost) where T

    #
    min_P = minimum(P_cost)
    max_P = maximum(P_cost)

    N_compounds = length(ΩS0_ppm)
    active_flags = Vector{BitVector}(undef, N_compounds)

    for n = 1:N_compounds

        N_systems = length(ΩS0_ppm[n])
        active_flags[n] = falses(N_systems)

        for i = 1:N_systems

            min_F = minimum(ΩS0_ppm[n][i])
            max_F = maximum(ΩS0_ppm[n][i])
            active_flags[n][i] = isoverlapped(min_F, max_F, min_P, max_P)
        end
    end

    return active_flags
end



function combineregionsresultstype1(As, hz2ppmfunc, P_y, cost_inds_set,
    w_ar,
    κs_β_an_ar::Vector{Vector{Vector{Vector{T}}}},
    d_an_ar::Vector{Vector{Vector{T}}},
    β_singlets_an_ar::Vector{Vector{Vector{T}}},
    d_singlets_an_ar::Vector{Vector{Vector{T}}}) where T

    N_compounds = length(As)
    N_regions = length(cost_inds_set)

    ΩS0 = getΩS(As)
    ΩS0_ppm = getPs(ΩS0, hz2ppmfunc)
    active_flags_ar = collect( findcoveredspinsystems(ΩS0_ppm, P_y[cost_inds_set[r]]) for r = 1:length(cost_inds_set) )

    κs_β_an = Vector{Vector{Vector{T}}}(undef, N_compounds)
    d_an = Vector{Vector{T}}(undef, N_compounds)
    β_singlets_an = Vector{Vector{T}}(undef, N_compounds)
    d_singlets_an = Vector{Vector{T}}(undef, N_compounds)
    w = Vector{T}(undef, N_compounds)

    for n = 1:N_compounds

        N_sys = length(As[n].αs)
        N_singlets = length(As[n].αs_singlets)

        # non-singlet spin systems.
        κs_β_an[n] = Vector{Vector{T}}(undef, N_sys)
        d_an[n] = Vector{T}(undef, N_sys)
        for i = 1:N_sys

            κs_β_an[n][i] = averageregions(κs_β_an_ar, n, i, i, active_flags_ar)
            d_an[n][i] = averageregions(d_an_ar, n, i, i, active_flags_ar)
        end

        β_singlets_an[n] = Vector{T}(undef, N_singlets)
        d_singlets_an[n] = Vector{T}(undef, N_singlets)
        for i = 1:N_singlets
            β_singlets_an[n][i] = averageregions(β_singlets_an_ar, n, i, i+N_sys, active_flags_ar)
            d_singlets_an[n][i] = averageregions(d_singlets_an_ar, n, i, i+N_sys, active_flags_ar)
        end

        N = sum( 1 for r = 1:N_regions if any(active_flags_ar[r][n]) )
        w[n] = sum( w_ar[r][n] ./ N for r = 1:N_regions if any(active_flags_ar[r][n]) )
    end

    return κs_β_an, d_an, β_singlets_an, d_singlets_an, w
end

# X[r][n][i][blah]. i is non-singlet spin system index.
# active_flags_ar[r][n][k]. k is all spin system index (including singlets).
# r is region index.
function averageregions(X::Vector{Vector{Vector{S}}},
    n::Int, i::Int, k::Int, active_flags_ar::Vector{Vector{BitVector}})::S where S

    N_regions = length(active_flags_ar)
    @assert length(X) == N_regions

    N = sum( 1 for r = 1:N_regions if active_flags_ar[r][n][k] )
    return sum( X[r][n][i] ./ N for r = 1:N_regions if active_flags_ar[r][n][k] )
end
