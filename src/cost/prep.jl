"""
region_min_dist is the minimum horizontal distance between regions, in ppm.
"""
function prepareoptim(config_path::String,
    molecule_names,
    hz2ppmfunc,
    U_cost0, y_cost0::Vector{Complex{T}},
    As;
    region_min_dist = 0.1) where T <: Real

    ## parse from config file.
    config_dict = Dict()
    if ispath(config_path)

        # TODO add error-handling if name is not found in the dictionary or filename does not exist.
        config_dict = JSON.parsefile(config_path)
    end

    N_compounds = length(molecule_names)

    λ_lbs = Vector{Vector{T}}(undef, N_compounds)
    λ_ubs = Vector{Vector{T}}(undef, N_compounds)
    Δsys_cs = Vector{Vector{T}}(undef, N_compounds)
    κs_β_orderings = Vector{Vector{Vector{Int}}}(undef, N_compounds)
    κs_β_DOFs = Vector{Vector{Int}}(undef, N_compounds)

    for n = 1:N_compounds

        dict = config_dict[molecule_names[n]] # TODO graceful error-handle.

        λ_lbs[n] = convert(Vector{T}, dict["λ_lb"])
        λ_ubs[n] = convert(Vector{T}, dict["λ_ub"])
        Δsys_cs[n] = convert(Vector{T}, dict["maximum chemical shift"])
        κs_β_orderings[n] = convert(Vector{Vector{Int}}, dict["κs_β ordering"])
        κs_β_DOFs[n] = convert(Vector{Int}, dict["κs_β degrees of freedom"])

    end

    # cs_delta_group = NMRSpecifyRegions.extractinfofromconfig( cs_config_path, molecule_names)
    # Δsys_cs = NMRSpecifyRegions.condenseΔcsconfig(cs_delta_group)

    ## get regions.
    ΩS0 = getΩS(As)
    ΩS0_ppm = getPs(ΩS0, hz2ppmfunc)

    exp_info = NMRSpecifyRegions.setupexperimentresults(molecule_names, ΩS0_ppm, Δsys_cs;
    min_dist = region_min_dist)

    P_cost0 = hz2ppmfunc.(U_cost0)
    cost_inds, cost_inds_set = NMRSpecifyRegions.getcostinds(exp_info, P_cost0)

    U_cost = U_cost0[cost_inds]
    P_cost = P_cost0[cost_inds]
    y_cost = y_cost0[cost_inds]

    return Δsys_cs, y_cost, U_cost, P_cost, exp_info, cost_inds, cost_inds_set,
    λ_lbs,
    λ_ubs,
    κs_β_orderings,
    κs_β_DOFs
end

"""
setupitpab(window::T, N_itp_samples::Int, domain_percentage::T;
   N_fit_positions::Int = 15,
   p0 = [0.5; 0.0],
   p_lb = [0.1; -5.0],
   p_ub = [0.6; 5.0],
   max_iters = 5000,
   xtol_rel = 1e-5,
   ftol_rel = 1e-5,
   maxtime = Inf,
   optim_algorithm = :LN_BOBYQA) where T <: Real

window ∈ (0,1)
optim_algorithm can be :GN_ESCH, :GN_ISRES, :LN_BOBYQA, :GN_DIRECT_L
"""
function setupitpab(window::T, N_itp_samples::Int, domain_percentage::T;
    N_fit_positions::Int = 15,
    p0 = [0.5; 0.0],
    p_lb = [0.1; -5.0],
    p_ub = [0.6; 5.0],
    max_iters = 5000,
    xtol_rel = 1e-5,
    ftol_rel = 1e-5,
    maxtime = Inf,
    optim_algorithm = :LN_BOBYQA) where T <: Real

    # get piece-wise linear monotone maps.
    infos, zs, p_range = MonotoneMaps.getendomorphismpiecewiselinear(zero(T), one(T),
    window;
    N_itp_samples = N_itp_samples,
    domain_percentage = domain_percentage)

    # get compact sigmoid parameters fitted to each of the piece-wise linear maps.
    # costfuncs, minxs, rets = MonotoneMaps.getcompactsigmoidparameters(infos;
    # p0 = p0, p_lb = p_lb, p_ub = p_ub, optim_algorithm = optim_algorithm,
    # max_iters = max_iters, N_fit_positions = N_fit_positions, xtol_rel = xtol_rel,
    # ftol_rel = ftol_rel, maxtime = maxtime)
    costfuncs, minxs, rets = MonotoneMaps.getcompactsigmoidparameters(infos;
    N_fit_positions = N_fit_positions, max_iters = max_iters,
    xtol_rel = xtol_rel, ftol_rel = ftol_rel, maxtime = maxtime,
    p0 = p0, p_lb = p_lb, p_ub = p_ub, optim_algorithm = optim_algorithm)
    #qs = collect( tt->MonotoneMaps.evalcompositelogisticprobit(tt, minxs[i][1], minxs[i][2]) for i = 1:length(minxs) )

    Δp = p_range[2]-p_range[1]
    itp_range = p_range[1]:Δp:p_range[end]

    # itp.
    a_samples = collect( minxs[i][1] for i = 1:length(minxs) )
    b_samples = collect( minxs[i][2] for i = 1:length(minxs) )

    a_itp = Interpolations.interpolate(a_samples, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    a_sitp = Interpolations.scale(a_itp, itp_range)
    a_setp = Interpolations.extrapolate(a_sitp, Interpolations.Flat()) # zero outside interp range.

    b_itp = Interpolations.interpolate(b_samples, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    b_sitp = Interpolations.scale(b_itp, itp_range)
    b_setp = Interpolations.extrapolate(b_sitp, Interpolations.Flat()) # zero outside interp range.

    return a_setp, b_setp,
        minxs, rets # debug
end
