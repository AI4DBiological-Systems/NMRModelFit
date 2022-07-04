

function loadregionsresultstype1(project_folder::String, filename::String)

    load_path = joinpath(project_folder, filename)
    dic = BSON.load(load_path)

    minf_ar = convert(Vector{Float64}, dic[:minf_ar])
    minx_ar = convert(Vector{Vector{Float64}}, dic[:minx_ar])
    ret_ar = convert(Vector{Symbol}, dic[:ret_ar])
    w_ar = convert(Vector{Vector{Float64}}, dic[:w_ar])
    κs_β_an_ar = convert(Vector{Vector{Vector{Vector{Float64}}}}, dic[:κs_β_an_ar])
    d_an_ar = convert(Vector{Vector{Vector{Float64}}}, dic[:d_an_ar])
    β_singlets_an_ar = convert(Vector{Vector{Vector{Float64}}}, dic[:β_singlets_an_ar])
    d_singlets_an_ar = convert(Vector{Vector{Vector{Float64}}}, dic[:d_singlets_an_ar])

    return minf_ar, minx_ar, ret_ar, w_ar, κs_β_an_ar, d_an_ar,
        β_singlets_an_ar, d_singlets_an_ar, dic
end

# mutates Bs.
function saveregionsresults!(Bs,
    minf_ar, minx_ar, ret_ar, w_ar::Vector{Vector{T}}, obj_func_ar,
    project_folder::String, filename::String) where T

    κs_β_an_ar, d_an_ar, β_singlets_an_ar,
    d_singlets_an_ar = NMRModelFit.exportβdfromregionsresults!(Bs, minx_ar, obj_func_ar)

    save_path = joinpath(project_folder, filename)
    BSON.bson(save_path,
    minf_ar = minf_ar,
    minx_ar = minx_ar,
    ret_ar = ret_ar,
    w_ar = w_ar,
    κs_β_an_ar = κs_β_an_ar,
    d_an_ar = d_an_ar,
    β_singlets_an_ar = β_singlets_an_ar,
    d_singlets_an_ar = d_singlets_an_ar)

    return nothing
end

function plotregion(P, U, q_U, P_y, y, P_cost, y_cost, display_threshold_factor, display_reduction_factor,
    save_folder, title_string, file_name;
    save_plot_flag = true,
    display_plot_flag = true,
    canvas_size = (1000, 400),
    spectrum_processing_func = real,
    ref_P,
    ref_y,
    ref_label_string = "reference")

    vertical_tag = ""
    if spectrum_processing_func == real
        vertical_tag = "real part"
    elseif spectrum_processing_func == imag
        vertical_tag = "imaginary part"
    elseif spectrum_processing_func == abs
        vertical_tag = "magnitude"
    end

    # # reduce the plotting positions for low signal regions. Otherwise the plot store size will be too large, and the time to load the plot will be long.
    # inds, _ = NMRSignalSimulator.prunelowsignalentries(q_U, display_threshold_factor, display_reduction_factor)
    # P_display = P[inds]
    # U_display = U[inds]
    # q_U_display = q_U[inds]
    P_display = P
    U_display = U
    q_U_display = q_U

    # plot.
    isdir(save_folder) || mkpath(save_folder)

    plots_save_path = joinpath(save_folder, file_name)
    #title_string = "$(project_name) alignment results, region $(r), real part"

    plot_obj = Plots.plot( P_display,
        spectrum_processing_func.(q_U_display),
        title = title_string,
        label = "model",
        seriestype = :line,
        ticks = :native,
        xlims = (P_display[1],P_display[end]),
        hover = P_display,
        linewidth = 4,
        xlabel = "ppm",
        ylabel = "$(vertical_tag) spectrum",
        xflip = true,
        size = canvas_size)

    Plots.plot!(plot_obj, P_y, spectrum_processing_func.(y), label = "full data",
        seriestype = :line,
        linestyle = :dot,
        xflip = true,
        linewidth = 4)

    Plots.plot!(plot_obj, P_cost, spectrum_processing_func.(y_cost), label = "fit data",
        markershape = :circle,
        seriestype = :scatter,
        xflip = true)

    if !isempty(ref_P) && !isempty(ref_y)
        Plots.plot!(plot_obj, ref_P, spectrum_processing_func.(ref_y), label = ref_label_string,
            seriestype = :line,
            linestyle = :dot,
            xflip = true,
            linewidth = 4)
    end

    if save_plot_flag
        Plots.savefig(plot_obj, plots_save_path)
    end

    if display_plot_flag
        display(plot_obj)
    end


    return nothing
end

function plotquantificationresults(As, Bs, ws, save_folder,
    P_y, y,
    P, U,
    obj_funcs, minfs, minxs, rets,
    display_reduction_factor, display_threshold_factor,
    cost_inds_set;
    canvas_size = (1000, 400),
    display_flag = false,
    save_plot_flag = true,
    N_viz = 50000,
    filename_prefix = "")

    U_rad = U .* (2*π)

    for r = 1:length(cost_inds_set)

        y_cost = y[cost_inds_set[r]]
        P_cost = P_y[cost_inds_set[r]]

        q2 = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs; w = ws[r])

        cost = obj_funcs[r](minxs[r])
        #ws[r][1] = 0.0
        q_U = q2.(U_rad)

        file_name = "$(filename_prefix)real_$(r).html"
        title_string = "$(project_name) region $(r), real part, cost = $(cost)"

        plotregion(P, U, q_U, P_y, y, P_cost, y_cost, display_threshold_factor, display_reduction_factor,
            save_folder, title_string, file_name;
            save_plot_flag = save_plot_flag,
            display_plot_flag = display_flag,
            canvas_size = canvas_size,
            spectrum_processing_func = real)

        #
        file_name = "$(filename_prefix)imaginary_$(r).html"
        title_string = "$(project_name) region $(r), imaginary part, cost = $(cost)"

        plotregion(P, U, q_U, P_y, y, P_cost, y_cost, display_threshold_factor, display_reduction_factor,
            save_folder, title_string, file_name;
            save_plot_flag = save_plot_flag,
            display_plot_flag = display_flag,
            canvas_size = canvas_size,
            spectrum_processing_func = imag)
    end

    return nothing
end



####
# taken from prep_full.jl
function setupfitsession(SH_config_path,
    surrogate_config_path,
    fit_config_path,
    H_params_path,
    dict_compound_to_filename,
    experiment_full_path,
    project_name,
    project_folder,
    molecule_names;
    u_offset = 0.5, # in ppm.
    solvent_ppm_guess = 4.7,
    solvent_window_ppm = 0.1,
    offset_ppm = 0.3,
    Δcs_max_scalar_default = 0.2,
    unique_cs_atol = 1e-6,
    prune_combo_Δc_bar_flag = true,
    region_min_dist = 0.1)


    ##### load experiment, normalize data.
    s_t, S, hz2ppmfunc, ppm2hzfunc, ν_0ppm, fs, SW, α_0ppm, β_0ppm, λ_0ppm, Ω_0ppm,
        results_0ppm, dic, α_solvent, β_solvent, λ_solvent, Ω_solvent,
        results_solvent = NMRDataSetup.loadspectrum(experiment_full_path;
        solvent_ppm = solvent_ppm_guess,
        solvent_window_ppm = solvent_window_ppm)
    #

    # start the first entry at the frequency corresponding to offset_Hz, which we set to 0.3 ppm.

    offset_Hz = ν_0ppm - (ppm2hzfunc(offset_ppm)-ppm2hzfunc(0.0))
    DFT_s = fft(s_t)
    U_DFT, U_y, U_inds = NMRDataSetup.getwraparoundDFTfreqs(length(s_t), fs, offset_Hz)

    S_U = DFT_s[U_inds]
    P_y = hz2ppmfunc.(U_y)

    val, ind = findmin( abs.(U_y .- ν_0ppm) )
    Z = abs(S_U[ind])

    y = S_U ./ Z
    #####
    #@assert 1==2

    ####### compute surrogate.
    λ0 = λ_0ppm

    # get a surrogate where K_{n,i} is encouraged to be no larger than `early_exit_part_size`.
    println("Timing: getphysicalparameters")
    @time Phys = NMRHamiltonian.getphysicalparameters(molecule_names,
        H_params_path,
        dict_compound_to_filename;
        unique_cs_atol = unique_cs_atol)

    println("Timing: setupmixtureSH()")
    @time As = NMRHamiltonian.setupmixtureSH(molecule_names,
        fs, SW, ν_0ppm,
        Phys;
        config_path = SH_config_path,
        prune_combo_Δc_bar_flag = prune_combo_Δc_bar_flag)


    ΩS_ppm = NMRModelFit.getΩSppm(As, hz2ppmfunc)
    ΩS_ppm_sorted = sort(NMRModelFit.combinevectors(ΩS_ppm))


    u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
    u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)



    println("fitproxies!()")
    Bs = NMRSignalSimulator.fitproxies(As, NMRSignalSimulator.SpinSysParamsType1(0.0), λ0;
        names = molecule_names,
        config_path = surrogate_config_path,
        Δcs_max_scalar_default = Δcs_max_scalar_default,
        u_min = u_min,
        u_max = u_max)

    #
    Bs2 = NMRSignalSimulator.fitproxies(As, NMRSignalSimulator.SpinSysParamsType2(0.0), λ0;
        names = molecule_names,
        config_path = surrogate_config_path,
        Δcs_max_scalar_default = Δcs_max_scalar_default,
        u_min = u_min,
        u_max = u_max)

    ####### end mixture proxy.


    ##### prepare fit positions.
    ### prepare positions.
    hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
    ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)
    P_y = hz2ppmfunc.(U_y)

    ΩS0 = NMRModelFit.getΩS(As)
    ΩS0_ppm = NMRModelFit.getPs(ΩS0, hz2ppmfunc)

    Δsys_cs, y_cost_all, U_cost_all, P_cost_all, exp_info, cost_inds, cost_inds_set,
        λ_lbs, λ_ubs, κs_β_orderings,
        κs_β_DOFs = NMRModelFit.prepareoptim(fit_config_path, molecule_names,
        hz2ppmfunc, U_y, y, As; region_min_dist = region_min_dist)
    #
    a_setp, b_setp, minxs,
        rets = NMRModelFit.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)

    return y, U_y, P_y, As, Bs, Bs2, fs, SW, ν_0ppm, cost_inds_set,
        Δsys_cs, u_min, u_max, cost_inds, Phys
end

function runfitregions(y::Vector{Complex{T}},
    U_y,
    P_y,
    P, U,
    As,
    Bs,
    fs,
    SW,
    Δsys_cs, a_setp, b_setp, cost_inds_set,
    u_min, u_max,
    Δcs_offset,
    Δcs_offset_singlets;
    maxeval = 50,
    save_BSON_flag = true,
    save_plot_flag = true,
    display_flag = true) where T <: Real

    if isempty(Δcs_offset) || isempty(Δcs_offset_singlets)
        Δcs_offset, Δcs_offset_singlets = NMRModelFit.initializeΔcstype1(Bs)
    end

    N_d = sum( NMRModelFit.getNdvars(Bs[n]) for n = 1:length(Bs) )

    shift_lb = -ones(N_d)
    shift_ub = ones(N_d)

    type_num = 1
    if typeof(Bs[1].ss_params) <: NMRModelFit.SpinSysParamsType2
        type_num = 2
    end

    println("Running fitregions()")
    obj_funcs, minfs, minxs, rets, ws = NMRModelFit.fitregions(y,
        U_y,
        P_y,
        As,
        Bs,
        fs,
        SW,
        Δsys_cs,
        a_setp, b_setp,
        shift_lb,
        shift_ub,
        cost_inds_set,
        Δcs_offset, Δcs_offset_singlets;
        N_starts = 100,
        local_optim_algorithm = NLopt.LN_BOBYQA,
        xtol_rel = 1e-9,
        maxeval = maxeval, # 2, # 50,
        maxtime = Inf,
        β_optim_algorithm = :GN_DIRECT_L,
        w_lb_default = 1e-1,
        w_ub_default = 100.0,
        β_max_iters = 500, # 2, # 500,
        β_xtol_rel = 1e-9,
        β_ftol_rel = 1e-9,
        β_maxtime = Inf);

    if save_BSON_flag
        save_name = "results_regions_type$(type_num).bson"
        saveregionsresults!(Bs, minfs, minxs, rets, ws,
        obj_funcs, project_folder, save_name)
    end

    display_reduction_factor = 1
    display_threshold_factor = 0.001/10
    if "L-Isoleucine" in molecule_names
        display_reduction_factor = 1
        display_threshold_factor = 0.001/10
    end

    plotquantificationresults(As, Bs, ws, project_folder,
        P_y, y, P, U,
        obj_funcs, minfs, minxs, rets,
        display_reduction_factor, display_threshold_factor,
        cost_inds_set;
        canvas_size = (1000, 400),
        display_flag = display_flag,
        save_plot_flag = save_plot_flag,
        N_viz = 50000,
        filename_prefix = "fit_type$(type_num)_")

    return obj_funcs, minfs, minxs, rets, ws
end

function preparetype2Δcsoffset(As::Vector{NMRModelFit.SHType{T}}, Phys,
    fs::T, SW::T, ν_0ppm,
    P_y,
    cost_inds_set,
    project_folder,
    type1_fit_result_filename) where T <: Real

    minf_ar, minx_ar, ret_ar, w_ar, κs_β_an_ar, d_an_ar,
        β_singlets_an_ar, d_singlets_an_ar,
        dic = loadregionsresultstype1(project_folder, type1_fit_result_filename)

    #
    cs_an, cs_singlets_an = NMRModelFit.getcs(As, Phys)
    Δcs_an_ravg, Δcs_an_singlets_ravg, Δcs_an_r, Δcs_singlets_an_r = NMRModelFit.getΔcstype2(P_y, cost_inds_set,
        cs_an, cs_singlets_an, d_an_ar, d_singlets_an_ar, fs, SW)


    return Δcs_an_ravg, Δcs_an_singlets_ravg, Δcs_an_r, Δcs_singlets_an_r
end

function runfitmixture(y_cost::Vector{Complex{T}}, P_cost, U_cost,
    As, Bs, fs, SW, Δsys_cs, Δcs_offset, Δcs_offset_singlets,
    a_setp, b_setp, P, U,
    project_name, save_folder, filename_prefix, BSON_filename;
    maxeval = 50,
    N_starts = 100,
    β_max_iters = 500,
    save_BSON_flag = true,
    save_plot_flag = true,
    display_plot_flag = true,
    N_viz = 50000,
    display_reduction_factor = 1,
    display_threshold_factor = 0.001/10,
    canvas_size = (1000, 400)) where T <: Real


    N_d = sum( NMRModelFit.getNdvars(Bs[n]) for n = 1:length(Bs) )

    shift_lb = -ones(N_d)
    shift_ub = ones(N_d)

    println("fitmodel")
    @time obj_func, minf, minx, ret,
    w = NMRModelFit.fitmodel(y_cost,
        U_cost,
        P_cost,
        As,
        Bs,
        fs,
        SW,
        Δsys_cs,
        a_setp, b_setp,
        shift_lb,
        shift_ub,
        Δcs_offset,
        Δcs_offset_singlets;
        N_starts = N_starts,
        local_optim_algorithm = NLopt.LN_BOBYQA,
        xtol_rel = 1e-9,
        maxeval = maxeval, # 2, # 50,
        maxtime = Inf,
        β_optim_algorithm = :GN_DIRECT_L,
        w_lb_default = 1e-1,
        w_ub_default = 100.0,
        β_max_iters = β_max_iters, # 2, # 500,
        β_xtol_rel = 1e-9,
        β_ftol_rel = 1e-9,
        β_maxtime = Inf)

    # @time obj_func, minf, minx, ret,
    # w = NMRModelFit.fitmodelBlackBoxOptim(y_cost,
    #     U_cost,
    #     P_cost,
    #     As,
    #     Bs,
    #     fs,
    #     SW,
    #     Δsys_cs,
    #     a_setp, b_setp,
    #     shift_lb,
    #     shift_ub,
    #     Δcs_offset,
    #     Δcs_offset_singlets;
    #     N_starts = N_starts,
    #     local_optim_algorithm = :adaptive_de_rand_1_bin,
    #     xtol_rel = 1e-9,
    #     maxeval = maxeval, # 2, # 50,
    #     maxtime = Inf,
    #     β_optim_algorithm = :GN_DIRECT_L,
    #     w_lb_default = 1e-1,
    #     w_ub_default = 100.0,
    #     β_max_iters = β_max_iters, # 2, # 500,
    #     β_xtol_rel = 1e-9,
    #     β_ftol_rel = 1e-9,
    #     β_maxtime = Inf)


    cost_obj_func = obj_func(minx)

    # load these cs into SH simulation to get proper amplitude.
    extracted_Δcs2 = ones(N_d) .* Inf
    NMRModelFit.extractmixtured!(extracted_Δcs2, Bs, As, 1, fs, SW)

    q2 = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs; w = w)

    U_cost_rad = U_cost .* (2*π)
    cost = sum( abs2.(q2.(U_cost_rad) - y_cost) ) # 1.6576135224225783 for serine-mod.
    println("abs(cost-cost_obj_func) = ", abs(cost-cost_obj_func))

    type_num = 1
    if typeof(Bs[1].ss_params) <: NMRModelFit.SpinSysParamsType2
        type_num = 2
    end

    if save_BSON_flag
        save_name = "mixture_fit_results_type$(type_num).bson"

        κs_β_an, d_an, κs_d_an, β_singlets_an, d_singlets_an = NMRModelFit.exportβdfromregionsresults!(Bs, minx, obj_func)

        save_path = joinpath(project_folder, BSON_filename)
        BSON.bson(save_path,
        minf = minf,
        minx = minx,
        ret = ret,
        w = w,
        κs_β_an = κs_β_an,
        κs_d_an = κs_d_an,
        d_an = d_an,
        β_singlets_an = β_singlets_an,
        d_singlets_an = d_singlets_an)

    end

    ### eval.
    U_rad = U .* (2*π)
    q_U = q2.(U_rad)

    #### initial.
    cost_initial = obj_func(zeros(T,length(minx)))
    println("initial cost = ", cost_initial)

    # load these cs into SH simulation to get proper amplitude.
    extracted_Δcs0 = ones(N_d) .* Inf
    NMRModelFit.extractmixtured!(extracted_Δcs0, Bs, As, 1, fs, SW)
    println("extracted_Δcs0 = ", extracted_Δcs0)

    q0 = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs; w = w)

    ### eval.
    q_U_initial = q0.(U_rad)

    ### plot.
    file_name = "$(filename_prefix)real.html"
    title_string = "$(project_name), cost = $(cost)"

    plotregion(P, U, q_U, P_y, y, P_cost, y_cost, display_threshold_factor, display_reduction_factor,
        save_folder, title_string, file_name;
        save_plot_flag = save_plot_flag,
        display_plot_flag = display_plot_flag,
        canvas_size = canvas_size,
        spectrum_processing_func = real,
        ref_P = P,
        ref_y = q_U_initial,
        ref_label_string = "initial model")
    #
    file_name = "$(filename_prefix)imaginary.html"
    title_string = "Fit for $(project_name), cost = $(cost)"

    plotregion(P, U, q_U, P_y, y, P_cost, y_cost, display_threshold_factor, display_reduction_factor,
        save_folder, title_string, file_name;
        save_plot_flag = save_plot_flag,
        display_plot_flag = display_plot_flag,
        canvas_size = canvas_size,
        spectrum_processing_func = imag,
        ref_P = P,
        ref_y = q_U_initial,
        ref_label_string = "initial model")

    return obj_func, minf, minx, ret, w, extracted_Δcs2
end
