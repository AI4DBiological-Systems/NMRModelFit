
function loadregion!(Bs::Vector{NMRModelFit.CompoundType{T, NMRModelFit.SpinSysParamsType1{T}}},
    project_folder::String) where T

    load_path = joinpath(project_folder, "region_$(r).bson")
    dic = BSON.load(load_path)
    minf = dic[:minf]
    minx = dic[:minx]
    ret = dic[:ret]
    w = dic[:w]

    for n = 1:length(Bs)

        if !isempty(dic[:κs_βs])

            tmp = convert(Vector{Vector{Float64}}, dic[:κs_βs][n])
            for i = 1:length(Bs[n].ss_params.κs_β)
                Bs[n].ss_params.κs_β[i][:] = tmp[i][:]
            end
            Bs[n].ss_params.d[:] = convert(Vector{Float64}, dic[:ds][n])

            if !isempty(Bs[n].β_singlets)
                Bs[n].β_singlets[:] = convert(Vector{Float64}, dic[:β_singletss][n])
                Bs[n].d_singlets[:] = convert(Vector{Float64}, dic[:d_singletss][n])
            end
        end
    end

    return minf, minx, rets, w
end

function saveregionresults(minfs, minxs, rets, ws,
    As, Bs, obj_funcs, project_folder::String)

    for r = 1:length(minxs)
        obj_funcs[r](minxs[r])

        save_path = joinpath(project_folder, "region_$(r).bson")
        BSON.bson(save_path, region_min_dist = region_min_dist,
        minf = minfs[r],
        minx = minxs[r],
        ret = rets[r],
        w = ws[r],
        κs_βs = collect( Bs[n].ss_params.κs_β for n = 1:length(Bs)),
        ds = collect( Bs[n].ss_params.d for n = 1:length(Bs)),
        β_singletss = collect( Bs[n].β_singlets for n = 1:length(Bs)),
        d_singletss = collect( Bs[n].d_singlets for n = 1:length(Bs)))
    end
end

function plotregion(P, U, q_U, P_y, y, P_cost, y_cost, display_threshold_factor, display_reduction_factor,
    save_folder, title_string, file_name;
    save_plot_flag = true,
    display_plot_flag = true,
    canvas_size = (1000, 400))

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
        real.(q_U_display),
        title = title_string,
        label = "model",
        seriestype = :line,
        ticks = :native,
        xlims = (P_display[1],P_display[end]),
        hover = P_display,
        linewidth = 4,
        xlabel = "ppm",
        ylabel = "real part of spectrum",
        xflip = true,
        size = canvas_size)

    Plots.plot!(plot_obj, P_y, real.(y), label = "full data",
        seriestype = :line,
        linestyle = :dot,
        xflip = true,
        linewidth = 4)

    Plots.plot!(plot_obj, P_cost, real.(y_cost), label = "fit data",
        markershape = :circle,
        seriestype = :scatter,
        xflip = true)

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
    region_min_dist,
    obj_funcs, minfs, minxs, rets,
    display_reduction_factor, display_threshold_factor,
    cost_inds_set, loop_range;
    canvas_size = (1000, 400),
    display_flag = false,
    save_plot_flag = true,
    N_viz = 50000)

    U = LinRange(u_min, u_max, N_viz)
    P = hz2ppmfunc.(U)
    U_rad = U .* (2*π)

    for r in loop_range

        y_cost = y[cost_inds_set[r]]
        P_cost = P_y[cost_inds_set[r]]

        q2 = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs; w = ws[r])

        cost = obj_funcs[r](minxs[r])
        #ws[r][1] = 0.0
        q_U = q2.(U_rad)

        file_name = "real_$(r).html"
        title_string = "$(project_name) fit results, region $(r), real part, cost = $(cost)"

        plotregion(P, U, q_U, P_y, y, P_cost, y_cost, display_threshold_factor, display_reduction_factor,
            save_folder, title_string, file_name;
            save_plot_flag = save_plot_flag,
            display_plot_flag = display_flag,
            canvas_size = canvas_size)
    end

    return nothing
end
