
function setupupdateLS(N_positions::Int, N_compounds::Int, y::Vector{Complex{T}}) where T <: Real

    F = Matrix{Complex{T}}(undef, N_positions, N_compounds)
    w = zeros(T, N_compounds)
    b = collect(reinterpret(T, y))

    return F, w, b
end

"""
F is an M x L complex-valued matrix, b is 2*M 1D array.
Updates solution to w (mutates it), and mutates intermediate variables F, b.
"""
function updatew!(  F::Matrix{Complex{T}},
    b,
    w::Vector{T},
    U_rad_LS,
    Bs,
    As,
    w_lower::Vector{T},
    w_upper::Vector{T}) where T <: Real

    @assert length(b) == 2*length(U_rad_LS)

    evaldesignmatrixw2!(F, U_rad_LS, Bs, As)

    ### LS solve.
    w[:], status_flag = solveBLScL(F, b, w_lower, w_upper)

    return status_flag
end

function solveBLScL( F::Matrix{Complex{T}},
    b,
    x_lower::Vector{T},
    x_upper::Vector{T}) where T <: Real

    # This is the formula.
    # BoundedLeastSquares.Quadratic(A'*A, A'*b)

    A = reinterpret(T, F)

    AtA = forcesymmetric(A'*A)
    Aty = A'*b

    Q = BoundedLeastSquares.Quadratic(AtA, Aty)

    x_opt, out_flag = BoundedLeastSquares.min_bound_constrained_quadratic(Q,
                        x_lower, x_upper)

    return x_opt, out_flag
end

function forcesymmetric(A::Matrix{T})::Matrix{T} where T <: Real
    return (A+A')./2
end

function evaldesignmatrixw2!(F::Matrix{Complex{T}},
    U_rad, Bs::Vector{CompoundType{T, SpinSysParamsType1{T}}}, As) where T <: Real

    M = length(U_rad)
    N = length(Bs)

    @assert size(F) == (M,N)

    fill!(F, NaN) # debug.

    for n = 1:N
        B = Bs[n]
        A = As[n]

        for m = 1:M

            # non-singlet spin systems.
            s = zero(Complex{T})
            for i in eachindex(B.qs)
                for k in eachindex(B.qs[i])
                    s += B.qs[i][k](U_rad[m] - B.ss_params.d[i], B.ss_params.κs_λ[i])
                end
            end
            F[m,n] = s

            # singlets.
            for k = 1:length(B.κs_λ_singlets)

                F[m,n] += NMRSpectraSimulator.evalsinglets(U_rad[m], B.d_singlets,
                A.αs_singlets, A.Ωs_singlets,
                B.β_singlets, B.λ0, B.κs_λ_singlets, B.κs_α_singlets)
            end

        end
    end

    return nothing
end

function evaldesignmatrixw2!(F::Matrix{Complex{T}},
    U_rad, Bs::Vector{CompoundType{T, SpinSysParamsType2{T}}}, As) where T <: Real

    M = length(U_rad)
    N = length(Bs)

    @assert size(F) == (M,N)

    fill!(F, NaN) # debug.

    for n = 1:N
        B = Bs[n]
        A = As[n]

        for m = 1:M

            # non-singlet spin systems.
            s = zero(Complex{T})
            for i in eachindex(B.qs)
                for k in eachindex(B.qs[i])
                    s += B.qs[i][k](U_rad[m] - B.ss_params.d[i][k], B.ss_params.κs_λ[i])
                end
            end
            F[m,n] = s

            # singlets.
            for k = 1:length(B.κs_λ_singlets)

                F[m,n] += NMRSpectraSimulator.evalsinglets(U_rad[m], B.d_singlets,
                A.αs_singlets, A.Ωs_singlets,
                B.β_singlets, B.λ0, B.κs_λ_singlets, B.κs_α_singlets)
            end

        end
    end

    return nothing
end
