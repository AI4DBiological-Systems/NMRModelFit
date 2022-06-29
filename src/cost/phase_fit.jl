function setupcostnesteddwarpw(Bs,
    As,
    fs::T,
    SW::T,
    LS_inds,
    U_rad_cost,
    y_cost::Vector{Complex{T}},
    Δsys_cs::Vector{Vector{T}},
    Δcs_offset::Vector{T},
    a_setp, b_setp;
    β_optim_algorithm = :GN_DIRECT_L,
    w_lb_default = 0.2,
    w_ub_default = 5.0,
    β_max_iters = 500,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf) where T <: Real


    ##### update functions.
    N_β = sum( getNβ(Bs[n]) for n = 1:length(Bs) )
    N_d = sum( getNd(Bs[n]) for n = 1:length(Bs) )
    p_buffer = zeros(T, N_d)

    st_ind_d = 1
    fin_ind_d = st_ind_d + N_d - 1
    updatedfunc = pp->updatemixturedwarp!(p_buffer, Bs, As, pp, st_ind_d, fs, SW,
        Δsys_cs, Δcs_offset, a_setp, b_setp)

    N_vars_set = [N_d; ]

    # β, κ update.
    run_optim, obj_func_β, E_BLS, w_BLS, b_BLS, updateβfunc, updatewfunc,
    q_β = setupβLSsolverw(β_optim_algorithm,
        Bs, As, LS_inds, U_rad_cost, y_cost;
        w_lb_default = w_lb_default,
        w_ub_default = w_ub_default,
        β_max_iters = β_max_iters,
        β_xtol_rel = β_xtol_rel,
        β_ftol_rel = β_ftol_rel,
        β_maxtime = β_maxtime)

    #### extract parameters from p.
    getshiftfunc = pp->pp[st_ind_d:fin_ind_d]

    # model. (optional)
    #f = uu->evalclproxymixture(uu, As, Bs; w = w_BLS)
    f = q_β

    return f, updatedfunc, getshiftfunc, N_vars_set,
    run_optim, obj_func_β, E_BLS, w_BLS, b_BLS, updateβfunc, updatewfunc, q_β
end


function setupβLSsolverw(optim_algorithm,
    Bs,
    As,
    LS_inds,
    U_rad_cost,
    y_cost::Vector{Complex{T}};
    w_lb_default = 0.2,
    w_ub_default = 5.0,
    β_max_iters = 50,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf) where T

    N_β = sum( getNβ(Bs[n]) for n = 1:length(Bs) )

    β_lb = ones(T, N_β) .* (-π)
    β_ub = ones(T, N_β) .* (π)

    p_lb = β_lb
    p_ub = β_ub

    q, updateβfunc, updatewfunc, E_BLS, w_BLS, b_BLS,
    getβfunc = setupcostβLSw(Bs, As, LS_inds, U_rad_cost, y_cost;
        w_lb_default = w_lb_default,
        w_ub_default = w_ub_default)

    # q is simulated spectra. f is cost function.
    f = pp->costβLSw(U_rad_cost, y_cost, updateβfunc, updatewfunc, pp,
    E_BLS, w_BLS, b_BLS, q)

    df = xx->FiniteDiff.finite_difference_gradient(f, xx)

    opt = NLopt.Opt(optim_algorithm, N_β)

    run_optim = pp->runNLopt!(opt,
        pp,
        f,
        df,
        p_lb,
        p_ub;
        max_iters = β_max_iters,
        xtol_rel = β_xtol_rel,
        ftol_rel = β_ftol_rel,
        maxtime = β_maxtime)

    return run_optim, f, E_BLS, w_BLS, b_BLS, updateβfunc, updatewfunc, q
end

function costβLSw(U_rad,
    S_U::Vector{Complex{T}},
    updateβfunc, updatewfunc,
    p::Vector{T},
    E_BLS::Matrix{Complex{T}}, w_BLS::Vector{T}, b_BLS,
    q)::T where T <: Real

    updateβfunc(p)
    updatewfunc(1.0)

    # ## l-2 costfunc.
    # cost = zero(T)
    # for m = 1:length(S_U)
    #
    #     q_u = q(U_rad[m])
    #
    #     cost += abs2( q_u - S_U[m] )
    # end
    #
    # # faster l-2 cost compute.
    # B = reinterpret(T, E_BLS)
    # tmp = B*w_BLS - b_BLS
    # cost = dot(tmp, tmp)
    cost = norm(reinterpret(T, E_BLS)*w_BLS - b_BLS)^2 # compact version.

    return cost
end


function setupcostβLSw(Bs::Vector{CompoundType{T, SST}},
    As,
    LS_inds,
    U0_rad,
    y0::Vector{Complex{T}};
    w_lb_default = 0.2,
    w_ub_default = 5.0) where {T <: Real, SST}

    N_β = sum( getNβ(Bs[n]) for n = 1:length(Bs) )

    st_ind_β = 1
    fin_ind_β = st_ind_β + N_β - 1
    updateβfunc = pp->updateβ!(Bs, pp, st_ind_β)

    ### LS w.
    U_rad_LS = U0_rad[LS_inds]

    N_w_vars = length(Bs)
    E_BLS, w_BLS, b_BLS = setupupdateLS(length(U_rad_LS), N_w_vars, y0[LS_inds])

    w_lb = ones(N_w_vars) .* w_lb_default
    w_ub = ones(N_w_vars) .* w_ub_default

    updatewfunc = xx->updatew!(E_BLS, b_BLS, w_BLS,
    U_rad_LS, Bs, As, w_lb, w_ub)

    #### extract parameters from p.
    getβfunc = pp->pp[st_ind_β:fin_ind_β]

    f = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs; w = w_BLS)

    return f, updateβfunc, updatewfunc, E_BLS, w_BLS, b_BLS, getβfunc
end

function updateβ!(Bs::Vector{CompoundType{T,SST}},
    p::Vector{T},
    st_ind::Int)::Int where {T <: Real, SST}

    j = st_ind - 1

    for n = 1:length(Bs)
        for i = 1:length(Bs[n].ss_params.κs_β)
            for l = 1:length(Bs[n].ss_params.κs_β[i])

                j += 1
                Bs[n].ss_params.κs_β[i][l] = p[j]
            end
        end

        for i = 1:length(Bs[n].β_singlets)
            j += 1
            Bs[n].β_singlets[i] = p[j]
        end
    end

    return j
end
