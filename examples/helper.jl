

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
    BSON.bson(save_path, region_min_dist = region_min_dist,
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
    N_viz = 50000,
    filename_prefix = "")

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

        file_name = "$(filename_prefix)real_$(r).html"
        title_string = "$(project_name) fit results, region $(r), real part, cost = $(cost)"

        plotregion(P, U, q_U, P_y, y, P_cost, y_cost, display_threshold_factor, display_reduction_factor,
            save_folder, title_string, file_name;
            save_plot_flag = save_plot_flag,
            display_plot_flag = display_flag,
            canvas_size = canvas_size)
    end

    return nothing
end
