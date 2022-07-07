function fitregions(y::Vector{Complex{T}}, U_y, P_y, As, Bs, fs, SW,
    Δsys_cs, a_setp, b_setp,
    shift_lb::Vector{T},
    shift_ub::Vector{T},
    cost_inds_set::Vector{Vector{Int}},
    Δcs_offset, Δcs_offset_singlets;
    N_starts = 100,
    local_optim_algorithm = NLopt.LN_BOBYQA,
    xtol_rel = 1e-3,
    maxeval = 50,
    maxtime = Inf,
    β_optim_algorithm = :GN_DIRECT_L,
    w_lb_default = 1e-2,
    w_ub_default = 1e2,
    β_max_iters = 500,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf) where T <: Real

    N_regions = length(cost_inds_set)
    minxs = Vector{Vector{T}}(undef, N_regions)
    rets = Vector{Symbol}(undef, N_regions)
    minfs = Vector{T}(undef, N_regions)
    obj_funcs = Vector{Function}(undef, N_regions)
    ws = Vector{Vector{T}}(undef, N_regions)

    for r = 1:length(cost_inds_set)
        println("Working on region $(r)")
        obj_funcs[r], minfs[r], minxs[r], rets[r],
        ws[r] = fitmodel(y[cost_inds_set[r]],
            U_y[cost_inds_set[r]],
            P_y[cost_inds_set[r]],
            As,
            Bs,
            fs,
            SW,
            Δsys_cs,
            a_setp, b_setp,
            shift_lb,
            shift_ub,
            Δcs_offset, Δcs_offset_singlets;
            N_starts = N_starts,
            local_optim_algorithm = local_optim_algorithm,
            xtol_rel = xtol_rel,
            maxeval = maxeval,
            maxtime = maxtime,
            β_optim_algorithm = β_optim_algorithm,
            w_lb_default = w_lb_default,
            w_ub_default = w_ub_default,
            β_max_iters = β_max_iters,
            β_xtol_rel = β_xtol_rel,
            β_ftol_rel = β_ftol_rel,
            β_maxtime = β_maxtime)
    end

    return obj_funcs, minfs, minxs, rets, ws
end

function fitmodel(y_cost::Vector{Complex{T}},
    U_cost,
    P_cost,
    As,
    Bs,
    fs,
    SW,
    Δsys_cs,
    a_setp, b_setp,
    shift_lb::Vector{T},
    shift_ub::Vector{T},
    Δcs_offset, Δcs_offset_singlets;
    LS_inds = 1:length(U_cost),
    N_starts = 100,
    local_optim_algorithm = NLopt.LN_BOBYQA,
    xtol_rel = 1e-3,
    maxeval = 50,
    maxtime = Inf,
    β_optim_algorithm = :GN_DIRECT_L,
    w_lb_default = 1e-3,
    w_ub_default = 1e2,
    β_max_iters = 500,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf) where T <: Real

    # prepare.
    N_d = sum( getNdvars(Bs[n]) for n = 1:length(Bs) )
    @assert length(shift_ub) == length(shift_lb) == N_d

    U_rad_cost = U_cost .* (2*π)

    # setup inner optim over β.
    q, updatedfunc, getshiftfunc, N_vars_set,
    run_optim, obj_func_β, E_BLS, w_BLS, b_BLS, updateβfunc, updatewfunc,
    q_β = setupcostnesteddwarpw(Bs, As, fs, SW, LS_inds, U_rad_cost,
        y_cost, Δsys_cs, Δcs_offset, Δcs_offset_singlets, a_setp, b_setp;
        β_optim_algorithm = β_optim_algorithm,
        w_lb_default = w_lb_default,
        w_ub_default = w_ub_default,
        β_max_iters = β_max_iters,
        β_xtol_rel = β_xtol_rel,
        β_ftol_rel = β_ftol_rel,
        β_maxtime = β_maxtime)

    # set up outer optim over shifts.
    N_β = sum( getNβ(Bs[n]) for n = 1:length(Bs) )

    p_β = zeros(T, N_β) # persistant buffer.

    obj_func = pp->costnesteddw(U_rad_cost, y_cost, updatedfunc, updatewfunc,
    pp, Bs, run_optim, E_BLS, w_BLS, b_BLS, p_β)

    # optim.
    prob = MultistartOptimization.MinimizationProblem(obj_func, shift_lb, shift_ub)

    local_method = MultistartOptimization.NLoptLocalMethod(; algorithm = local_optim_algorithm,
    xtol_rel = xtol_rel,
    maxeval = maxeval,
    maxtime = maxtime)

    println("N_starts = ", N_starts)

    multistart_method = MultistartOptimization.TikTak(N_starts)
    ret_mo = MultistartOptimization.multistart_minimization(multistart_method,
        local_method, prob)
    #
    println("ret_mo = ", ret_mo)
    minf = ret_mo.value
    minx = ret_mo.location
    ret = :None
    if hasproperty(ret_mo, :ret)
        ret = ret_mo.ret
    end

    # force w_BLS to update.
    obj_func(minx)

    return obj_func, minf, minx, ret, w_BLS
end

"""
Note that NLopt handle exceptions gracefully. Check the termination status to see if an exception occurs (would return FORCED_STOP)
"""
function costnesteddw(U,
    S_U::Vector{Complex{T}},
    updatedfunc,
    updatewfunc,
    p::Vector{T}, Bs,
    run_optim_β_κ::Function,
    E_BLS::Matrix{Complex{T}}, w_BLS::Vector{T}, b_BLS,
    p_β::Vector{T})::T where T <: Real

    updatedfunc(p)

    ### minimize inner problem.
    fill!(p_β, zero(T)) # always start from 0-phase?
    minf, minx, ret, N_evals = run_optim_β_κ(p_β)
    p_β[:] = minx # take out?

    # ensure Bs is updated with the latest β.
    #updateβ!(Bs, κs_β_orderings, κs_β_DOFs, minx, 1)
    updateβ!(Bs, minx, 1)
    updatewfunc(1.0)

    # evaluate cost.
    cost = norm(reinterpret(T, E_BLS)*w_BLS - b_BLS)^2

    return cost
end


# #### black box optim.
# function fitmodelBlackBoxOptim(y_cost::Vector{Complex{T}},
#     U_cost,
#     P_cost,
#     As,
#     Bs,
#     fs,
#     SW,
#     Δsys_cs,
#     a_setp, b_setp,
#     shift_lb::Vector{T},
#     shift_ub::Vector{T},
#     Δcs_offset, Δcs_offset_singlets;
#     initial_guess::Vector{T} = zeros(T, length(shift_lb)),
#     LS_inds = 1:length(U_cost),
#     N_starts = 100,
#     local_optim_algorithm = :adaptive_de_rand_1_bin,
#     xtol_rel = 1e-3,
#     maxeval = 50,
#     maxtime = Inf,
#     β_optim_algorithm = :GN_DIRECT_L,
#     w_lb_default = 1e-3,
#     w_ub_default = 1e2,
#     β_max_iters = 500,
#     β_xtol_rel = 1e-9,
#     β_ftol_rel = 1e-9,
#     β_maxtime = Inf) where T <: Real
#
#     # prepare.
#     N_d = sum( getNdvars(Bs[n]) for n = 1:length(Bs) )
#     @assert length(shift_ub) == length(shift_lb) == N_d
#
#     U_rad_cost = U_cost .* (2*π)
#
#     # setup inner optim over β.
#     q, updatedfunc, getshiftfunc, N_vars_set,
#     run_optim, obj_func_β, E_BLS, w_BLS, b_BLS, updateβfunc, updatewfunc,
#     q_β = setupcostnesteddwarpw(Bs, As, fs, SW, LS_inds, U_rad_cost,
#         y_cost, Δsys_cs, Δcs_offset, Δcs_offset_singlets, a_setp, b_setp;
#         β_optim_algorithm = β_optim_algorithm,
#         w_lb_default = w_lb_default,
#         w_ub_default = w_ub_default,
#         β_max_iters = β_max_iters,
#         β_xtol_rel = β_xtol_rel,
#         β_ftol_rel = β_ftol_rel,
#         β_maxtime = β_maxtime)
#
#     # set up outer optim over shifts.
#     N_β = sum( getNβ(Bs[n]) for n = 1:length(Bs) )
#
#     p_β = zeros(T, N_β) # persistant buffer.
#
#     obj_func = pp->costnesteddw(U_rad_cost, y_cost, updatedfunc, updatewfunc,
#     pp, Bs, run_optim, E_BLS, w_BLS, b_BLS, p_β)
#
#     # optim.
#     search_range = collect( (shift_lb[i], shift_ub[i]) for i = 1:length(shift_lb) )
#
#     res = BlackBoxOptim.bboptimize(obj_func, initial_guess;
#         SearchRange = search_range,
#         Method = local_optim_algorithm,
#         MaxTime = maxtime,
#         MaxFuncEvals = maxeval)
#
#     minx = BlackBoxOptim.best_candidate(res)
#     minf = BlackBoxOptim.best_fitness(res)
#     ret = Symbol(res.stop_reason)
#
#     # force w_BLS to update.
#     obj_func(minx)
#
#     return obj_func, minf, minx, ret, w_BLS
# end
