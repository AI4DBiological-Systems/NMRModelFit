


function getNd(A::CompoundType{T,SST})::Int where {T,SST}

    counter_sys = 0
    for i = 1:length(A.ss_params.d)
        counter_sys += length(A.ss_params.d[i])
    end

    return counter_sys + length(A.d_singlets)
end

function getNβ(A::CompoundType{T,SST}) where {T,SST}

    counter_sys = 0
    for i = 1:length(A.ss_params.κs_β)
        counter_sys += length(A.ss_params.κs_β[i])
    end

    return counter_sys + length(A.β_singlets)
end

"""
convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where T <: Real
returns (x-a)*(d-c)/(b-a)+c

converts x ∈ [a,b] to y ∈ [c,d].
"""
function convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where T <: Real

    return (x-a)*(d-c)/(b-a)+c
end

function getΩSppm(As::Vector{NMRHamiltonian.SHType{T}}, hz2ppmfunc) where T

    ΩS_ppm = Vector{Vector{T}}(undef, length(As))

    for (n,A) in enumerate(As)

        ΩS_ppm[n] = hz2ppmfunc.( combinevectors(A.Ωs) ./ (2*π) )

        tmp = hz2ppmfunc.( A.Ωs_singlets ./ (2*π) )
        push!(ΩS_ppm[n], tmp...)
    end

    return ΩS_ppm
end

function combinevectors(x::Vector{Vector{T}})::Vector{T} where T

    if isempty(x)
        return Vector{T}(undef, 0)
    end

    N = sum(length(x[i]) for i = 1:length(x))

    y = Vector{T}(undef,N)

    st_ind = 0
    fin_ind = 0
    for i = 1:length(x)
        st_ind = fin_ind + 1
        fin_ind = st_ind + length(x[i]) - 1

        y[st_ind:fin_ind] = x[i]
    end

    return y
end

"""
getΩS(As::Vector{CompoundFIDType{T,SST}}) where {T,SST}

returns ΩS::Vector{Vector{Vector{T}}}
"""
function getΩS(As::Vector{NMRHamiltonian.SHType{T}}) where T

    ΩS = Vector{Vector{Vector{T}}}(undef, length(As))

    for n = 1:length(As)

        ΩS[n] = Vector{Vector{T}}(undef, length(As[n].Ωs) + length(As[n].Ωs_singlets))
        for i = 1:length(As[n].Ωs)

            ΩS[n][i] = copy(As[n].Ωs[i])

        end

        for i = 1:length(As[n].Ωs_singlets)
            ΩS[n][length(As[n].Ωs)+i] = [ As[n].Ωs_singlets[i]; ]
        end
    end

    return ΩS
end

"""
getPs( ΩS::Vector{Vector{Vector{T}}}, hz2ppmfunc) where T <: Real

returns Ps::Vector{Vector{Vector{T}}}
"""
function getPs( ΩS::Vector{Vector{Vector{T}}}, hz2ppmfunc) where T <: Real

    N_compounds = length(ΩS)

    Ps = Vector{Vector{Vector{T}}}(undef, N_compounds)
    for n = 1:N_compounds

        Ps[n] = Vector{Vector{T}}(undef, length(ΩS[n]))
        for i = 1:length(ΩS[n])

            Ps[n][i] = Vector{T}(undef, length(ΩS[n][i]))
            for l = 1:length(ΩS[n][i])

                Ps[n][i][l] = hz2ppmfunc( ΩS[n][i][l]/(2*π) )
            end
        end
    end

    return Ps
end





function runNLopt!(  opt,
    p0::Vector{T},
    obj_func,
    grad_func,
    p_lbs,
    p_ubs;
    max_iters = 10000,
    xtol_rel = 1e-12,
    ftol_rel = 1e-12,
    maxtime = Inf) where T

    @assert length(p0) == length(p_lbs) == length(p_ubs)

    opt.maxeval = max_iters
    opt.lower_bounds = p_lbs
    opt.upper_bounds = p_ubs
    opt.xtol_rel = xtol_rel
    opt.ftol_rel = ftol_rel
    opt.maxtime = maxtime


    opt.min_objective = (xx, gg)->genericNLoptcostfunc!(xx, gg, obj_func, grad_func)

    # optimize.
    (minf, minx, ret) = NLopt.optimize(opt, p0)

    N_evals = opt.numevals

    return minf, minx, ret, N_evals
end

function genericNLoptcostfunc!(x::Vector{T}, df_x, f, df)::T where T <: Real

    #
    if length(df_x) > 0
        df_x[:] = df(x)
    end

    return f(x)
end
