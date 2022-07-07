
function evalclsinglets(u_rad::T, Ωs, λs, z::Vector{Complex{T}}) where T <: Real
    @assert length(Ωs) == length(λs) == length(z)
    
    out = zero(Complex{T})
    for l = 1:length(Ωs)

        out += z[l]/(λs[l]+im*(u_rad-Ωs[l]))
    end

    return out
end

function updatez!(A::Matrix{Complex{T}},
    z::Vector{Complex{T}},
    b::Vector{Complex{T}},
    U_rad,
    Ωs, λs) where T <: Real

    evaldesignmatrixcL!(A, U_rad, Ωs, λs)
    z[:] = A\b

    return nothing
end

function evaldesignmatrixcL!(A::Matrix{Complex{T}},
    U_rad, Ωs::Vector{T}, λs::Vector{T}) where T <: Real

    M = length(U_rad)
    N = length(Ωs)
    @assert length(λs) == N

    @assert size(A) == (M,N)

    fill!(A, NaN) # debug.

    for n = 1:N
        λ = λs[n]
        Ω = Ωs[n]

        for m = 1:M
            A[m,n] = 1/(λ+im*(U_rad[m]-Ω))
        end
    end

    return nothing
end

function fitsinglets(y::Vector{Complex{T}},
    U_rad,
    Ω_lb::Vector{T}, Ω_ub, λ_lb::Vector{T}, λ_ub::Vector{T},
    Ωs_initial, λs_initial;
    max_iters = 2000,
    xtol_rel = 1e-9,
    ftol_rel = 1e-9,
    maxtime = Inf,
    optim_algorithm = :GN_DIRECT_L) where T

    N_singlets = length(Ω_lb)
    A = Matrix{Complex{Float64}}(undef, length(U_rad), N_singlets)
    z = Vector{Complex{Float64}}(undef, N_singlets)
    updatezfunc = (pp,qq)->updatez!(A, z, y, U_rad, pp, qq)

    f = xx->evalcostsinglets(xx, y, A, z, updatezfunc)
    df = xx->FiniteDiff.finite_difference_gradient(f, xx)

    N_vars = 2*N_singlets
    opt = NLopt.Opt(optim_algorithm, N_vars)

    minf, minx, ret, N_evals = runNLopt!(opt,
        [Ωs_initial; λs_initial],
        f,
        df,
        [Ω_lb; λ_lb],
        [Ω_ub; λ_ub];
        max_iters = max_iters,
        xtol_rel = xtol_rel,
        ftol_rel = ftol_rel,
        maxtime = maxtime)

    return minf, minx, ret, N_evals, f, A, z
end

function evalcostsinglets( p::Vector{T},
    y::Vector{Complex{T}},
    A::Matrix{Complex{T}},
    z::Vector{Complex{T}},
    updatezfunc) where T <: Real

    # parse.
    L = div(length(p),2)
    Ωs = p[1:L]
    λs = p[L+1:end]

    updatezfunc(Ωs, λs)

    # evaluate cost.
    cost = norm(A*z - y)^2

    return cost
end
