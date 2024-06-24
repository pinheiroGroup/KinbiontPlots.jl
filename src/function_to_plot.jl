using Plots
using CSV

#plot a file



"""

    plot_data(
    label_exp::String, 
    path_to_data::String; 
    path_to_annotation::Any = missing,
    path_to_plot="NA",
    display_plots=true,
    save_plots=false,
    overlay_plots=true, 
    do_blank_subtraction="NO", vg)
    avg_replicate=false, 
    correct_negative="thr_correction", 
    thr_negative=0.01 ,
    blank_value = 0.0,
    blank_array = [0.0],
    text_size = 10,
    x_size =
    y_size =
    )

This function plot all the data from .csv file, note that assume that the first colums is the time

# Arguments:
- `label_exp::String`: The label of the experiment.  
- `path_to_data::String`: The path to the .csv of data

# Key Arguments:
-  `path_to_annotation::Any = missing`: The path to the .csv of annotation 
- `path_to_plot= "NA"`:String, path to save the plots.
-  `save_plot=false` :Bool, save the plot or not
- ` display_plots=true`:Bool,  Whether or not diplay the plot in julia
-   `overlay_plots =true` :Bool, if true it does one plot overlaying  all dataset curves
- `blank_subtraction="NO"`: String, how perform the blank subtration, options "NO","avg_subtraction" (subtration of average value of blanks) and "time_avg" (subtration of  time average value of blanks).  
- `average_replicate=false` Bool, perform or not the average of replicates. Works only if an annotation path is provided
- `blank_value = 0.0`: used only if `path_to_annotation = missing`and `blank_subtraction != "NO "`. It is used as average value of the blank
- `blank_array = [0.0]`:used only if `path_to_annotation = missing`and `blank_subtraction != "NO "`. It is used as array of the blanks values
-  `correct_negative="thr_correction"`  : String, How to treat negative values after blank subtraction. If `"thr_correction"` it put a thr on the minimum value of the data with blank subracted,
    if `"blank_correction"` uses blank distribution to impute negative values, if `"remove"` the values are just removed.
- `x_size=300` : X Size of the plot,
- `y_size=300` : y Size of the plot,
- `text_size = 10` : size of the text in the plots

# Output:
- For this function the output are saved or displayed depending on the values of key arguments.
"""
function plot_data(
    label_exp::String, #label of the experiment
    path_to_data::String; # path to the folder to analyze
    path_to_annotation::Any = missing,# path to the annotation of the wells
    path_to_plot="NA", # path where to save Plots
    display_plots=true,# display plots in julia or not
    save_plots=false, # save the plot or not
    overlay_plots=true, # true a single plot for all dataset false one plot per well
    do_blank_subtraction="NO", # string on how to use blank (NO,avg_subtraction,time_avg)
    avg_replicate=false, # if true the average between replicates
    correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01 ,
    blank_value = 0.0,
    blank_array = [0.0],
    x_size=300,
    y_size =300,
    guidefontsize=18,
    tickfontsize=16,
    legendfontsize=10
    )
  
    names_of_annotated_df,properties_of_annotation,list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)
    # reading files
    dfs_data = CSV.File(path_to_data)

    # shaping df for the inference
    names_of_cols = propertynames(dfs_data)

    # excluding blank data and discarded wells
    if length(list_of_blank) > 0
        names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
    end

    if length(list_of_discarded) > 0
        names_of_cols = filter!(e -> !(e in list_of_discarded), names_of_cols)
    end

    times_data = dfs_data[names_of_cols[1]]
    if length(list_of_blank) > 0
        blank_array = reduce(vcat, [(dfs_data[k]) for k in list_of_blank])
        blank_array = convert(Vector{Float64}, blank_array)

        blank_value = blank_subtraction(
            dfs_data,
            list_of_blank;
            method=do_blank_subtraction
        )

    end


    ## considering replicates
    list_replicate = unique(properties_of_annotation)
    list_replicate = filter!(e -> e != "b", list_replicate)

    if avg_replicate == true

        dfs_data, names_of_cols = average_replicate(dfs_data, times_data, properties_of_annotation, names_of_annotated_df)


    end
    # creating the folder of data if one wants to save
    if save_plots == true
        mkpath(path_to_plot)
    end

    for well_name in names_of_cols[2:end]
        name_well = string(well_name)

        if avg_replicate == true
            data_values = copy(dfs_data[!, well_name])
        else
            data_values = copy(dfs_data[well_name])
        end

        # blank subtraction
        data_values = data_values .- blank_value
        index_missing = findall(ismissing, data_values)
        index_tot =  eachindex(data_values)
        index_tot =  setdiff(index_tot,index_missing)
        data = Matrix(transpose(hcat(times_data[index_tot], data_values[index_tot])))
        # correcting negative values after blank subtraction
        data = negative_value_correction(data,
            blank_array;
            method=correct_negative,
            thr_negative=thr_negative,)

        if display_plots
            if_display = display
        else
            if_display = identity
        end

        if overlay_plots
            if well_name == names_of_cols[2]
                if_display(
                    Plots.plot(
                        data[1, :],
                        data[2, :],
                        xlabel="Time",
                        ylabel="Arb. Units",
                        size=(y_size,x_size),
                        label=[name_well],
                        title=string(label_exp),
                        legend=:outertopright,
                        guidefontsize=guidefontsize,
                        tickfontsize=tickfontsize,
                        legendfontsize=legendfontsize,
                    ),
                )
            else
                if_display(
                    Plots.plot!(
                        data[1, :],
                        data[2, :],
                        xlabel="Time",
                        ylabel="Arb. Units",
                        label=[name_well],
                        size=(y_size,x_size),
                        title=string(label_exp),
                        legend=:outertopright,
                        guidefontsize=guidefontsize,
                        tickfontsize=tickfontsize,
                        legendfontsize=legendfontsize,
                    ),
                )
            end

            if save_plots
                png(string(path_to_plot, label_exp, ".png"))
            end
        else
            # save & not plot single plot
            if_display(
                Plots.plot(
                    data[1, :],
                    data[2, :],
                    xlabel="Time",
                    ylabel="Arb. Units",
                    label=["Data " nothing],
                    color=:black,
                    title=string(label_exp, " ", name_well),
                    guidefontsize=guidefontsize,
                    tickfontsize=tickfontsize,
                    legendfontsize=legendfontsize,
                    size=(y_size,x_size),

                ),
            )
            if save_plots
                png(string(path_to_plot, label_exp, "_", name_well, ".png"))
            end
        end

    end
end








function plot_fit_of_file(
        Kimchi_results::Any,
        path_to_data::String; # path to the folder to analyze
        path_to_annotation::Any = missing,# path to the annotation of the wells
        path_to_plot="NA", # path where to save Plots
        display_plots=true,# display plots in julia or not
        save_plots=false, # save the plot or not
        do_blank_subtraction="avg_subtraction", # string on how to use blank (NO,avg_subtraction,time_avg)
        avg_replicate=false, # if true the average between replicates
        correct_negative="thr_correction", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
        thr_negative=0.01 ,
        blank_value = 0.0,
        blank_array = [0.0],
        pt_smoothing_derivative = 7,
        x_size=300,
        y_size =300,
        guidefontsize=18,
        tickfontsize=16,
        legendfontsize=10
        )
      


        Kimchi_results_matrix = Kimchi_results[2]
        Kimchi_method = Kimchi_results[1]



        names_of_annotated_df,properties_of_annotation,list_of_blank, list_of_discarded = reading_annotation(path_to_annotation)
        # reading files
        dfs_data = CSV.File(path_to_data)
    
        # shaping df for the inference
        names_of_cols = propertynames(dfs_data)
    
        # excluding blank data and discarded wells
        if length(list_of_blank) > 0
            names_of_cols = filter!(e -> !(e in list_of_blank), names_of_cols)
        end
    
        if length(list_of_discarded) > 0
            names_of_cols = filter!(e -> !(e in list_of_discarded), names_of_cols)
        end
    
        times_data = dfs_data[names_of_cols[1]]

        if length(list_of_blank) > 0
            blank_array = reduce(vcat, [(dfs_data[k]) for k in list_of_blank])
            blank_array = convert(Vector{Float64}, blank_array)
    
            blank_value = blank_subtraction(
                dfs_data,
                list_of_blank;
                method=do_blank_subtraction
            )
    
        end
    
    
        ## considering replicates
        list_replicate = unique(properties_of_annotation)
        list_replicate = filter!(e -> e != "b", list_replicate)
    
        if avg_replicate == true
    
            dfs_data, names_of_cols = average_replicate(dfs_data, times_data, properties_of_annotation, names_of_annotated_df)
    
    
        end
        # creating the folder of data if one wants to save
        if save_plots == true
            mkpath(path_to_plot)
        end
    
        for well_name in names_of_cols[2:end]
            name_well = string(well_name)


            if avg_replicate == true
                data_values = copy(dfs_data[!, well_name])
            else
                data_values = copy(dfs_data[well_name])
            end
    
            # blank subtraction
            data_values = data_values .- blank_value

            index_missing = findall(ismissing, data_values)
            index_tot =  eachindex(data_values)
            index_tot =  setdiff(index_tot,index_missing)
            data = Matrix(transpose(hcat(times_data[index_tot], data_values[index_tot])))
            # correcting negative values after blank subtraction
            data = negative_value_correction(data,
                blank_array;
                method=correct_negative,
                thr_negative=thr_negative,)
    
            if display_plots
                if_display = display
            else
                if_display = identity
            end
    
          index_of_well = findfirst(Kimchi_results_matrix[2,:] .== string(well_name) )

           # plotting fit 
          results_of_specific_well = Kimchi_results_matrix[:,index_of_well]

          if Kimchi_method == "Log-Lin"

            plot_log_lin(save_plots,
            if_display,
            data,
            results_of_specific_well;
            guidefontsize=guidefontsize,
            tickfontsize=tickfontsize,
            legendfontsize=legendfontsize,
            y_size =y_size,
            x_size =x_size,
            pt_smoothing_derivative = pt_smoothing_derivative )

          elseif Kimchi_method == "ODE"
            
            plot_ode_fit(save_plots,if_display,data,results_specific_well)

          elseif Kimchi_method == "ODE_segmentation"

            plot_ode_fit(save_plots,if_display,data,results_specific_well)


          elseif Kimchi_method == "segment_analysis"

            plot_seg_fit(save_plots,if_display,data,results_specific_well)


            


         elseif Kimchi_method == "NL"

            plot_NL_fit(save_plots,if_display,data,results_specific_well)


                
         elseif Kimchi_method == "NL_model_selection"

            plot_NL_fit(save_plots,if_display,data,results_specific_well)


         elseif Kimchi_method == "NL_segmentation"

            plot_NL_fit(save_plots,if_display,data,results_specific_well)



         end   






    
        end
end    


function plot_log_lin(save_plots,
    if_display,
    data,
    results_specific_well;
    guidefontsize=18,
    tickfontsize=12,
    legendfontsize=10,
    y_size = 300,
    x_size = 300,
    pt_smoothing_derivative = 7
    )

    name_well =   results_specific_well[2]
    label_exp =   results_specific_well[1]


    coeff_2 =  results_specific_well[7]
    coeff_1 = results_specific_well[12]

    index_of_t_start =  findfirst(data[1, :] .>=  results_specific_well[3])
    index_of_t_end =  findfirst(data[1, :] .>= results_specific_well[4])

    data_to_fit_times = data[1, index_of_t_start:index_of_t_end]
    data_to_fit_values = log.(data[2, index_of_t_start:index_of_t_end])

    N = length(data_to_fit_times)
    M = [ones(N) data_to_fit_times]
    mean_x = mean(data_to_fit_times)
    sigma_a = sigma_b = r = zeros(N)
    Theoretical_fitting = coeff_1 .+ data_to_fit_times .* coeff_2

    Cantrell_errors = sqrt(sum((data_to_fit_values - coeff_2 * data_to_fit_times .- coeff_1) .^ 2) / (N - 2))  # goodness of fit
    sigma_b = sqrt(1 / sum((data_to_fit_times .- mean_x) .^ 2))
    sigma_a = Cantrell_errors * sqrt(1 / N + mean_x^2 * sigma_b^2)
    sigma_b *= Cantrell_errors
    # Pearson's correlation coefficient
    d = TDist(N - 2)     # t-Student distribution with N-2 degrees of freedom
    cf = quantile(d, 0.975)  # correction factor for 95% confidence intervals (two-tailed distribution)
    confidence_band = cf * Cantrell_errors * sqrt.(1 / N .+ (data_to_fit_times .- mean(data_to_fit_times)) .^ 2 / var(data_to_fit_times) / (N - 1))

    if_display(
        Plots.scatter(
            data[1, :],
            log.(data[2, :]),
            xlabel="Time",
            ylabel="Log(Arb. Units)",
            label=["Data " nothing],
            markersize=1,
            color=:black,
            title=string(label_exp, " ", name_well),
            guidefontsize=guidefontsize,
            tickfontsize=tickfontsize,
            legendfontsize=legendfontsize,
            size=(y_size,x_size),

        ),
    )

    if_display(
        Plots.plot!(
            data_to_fit_times,
            Theoretical_fitting,
            ribbon=confidence_band,
            xlabel="Time ",
            ylabel="Log(Arb. Units)",
            label=[string("Fitting Log-Lin ") nothing],
            c=:red,
            guidefontsize=guidefontsize,
            tickfontsize=tickfontsize,
            legendfontsize=legendfontsize,
            size=(y_size,x_size),

        ),
    )
    if_display(
        Plots.vline!(
            [data_to_fit_times[1], data_to_fit_times[end]],
            c=:black,
            label=[string("Window of exp. phase ") nothing],
            guidefontsize=guidefontsize,
            tickfontsize=tickfontsize,
            legendfontsize=legendfontsize,
            size=(y_size,x_size),

        ),
    )
    if save_plots
        png(string(path_to_plot, label_exp, "_Log_Lin_Fit_", name_well, ".png"))
    end

    specific_gr = Kimchi.specific_gr_evaluation(data, pt_smoothing_derivative)

    specific_gr_times = [
        (data[1, r] + data[1, (r+pt_smoothing_derivative)]) / 2 for
        r = 1:1:(eachindex(data[2, :])[end].-pt_smoothing_derivative)
    ]

    if_display(
        Plots.scatter(
            specific_gr_times,
            specific_gr,
            xlabel="Time ",
            ylabel="1 /time ",
            label=[string("Dynamics growth rate ") nothing],
            c=:red,
            guidefontsize=guidefontsize,
            tickfontsize=tickfontsize,
            legendfontsize=legendfontsize,
            size=(y_size,x_size),

        ),
    )
    if_display(
        Plots.vline!(
            [data_to_fit_times[1], data_to_fit_times[end]],
            c=:black,
            label=[string("Window of exp. phase ") nothing],
               guidefontsize=guidefontsize,
        tickfontsize=tickfontsize,
        legendfontsize=legendfontsize,
        size=(y_size,x_size),
    ) )
    if save_plots
        png(string(path_to_plot, label_exp, "_specific_gr_dynamics_", name_well, ".png"))
    end

end

function plot_ode_fit(save_plots,
    if_display,
    data,
    results_specific_well;
    guidefontsize=18,
    tickfontsize=12,
    legendfontsize=10,
    y_size = 300,
    x_size = 300,
    pt_avg = 3,
    smoothing =false
    )
    name_well =   results_specific_well[2]
    label_exp =   results_specific_well[1]
    max_t = data[1, end]
    min_t = data[1, 1]
    tspan = (min_t, max_t)
    tsteps = data[1, :]

    model_string = results_specific_well[3]
    param_array =results_specific_well[3: (end - 4)]
    u0 = generating_IC(data, model, smoothing, pt_avg)
    ODE_prob = model_selector(model, u0, tspan)

    remade_solution = solve(remake(ODE_prob, p=res.u), integrator, saveat=tsteps)
    sol_time = reduce(hcat, remade_solution.t)
    sol_fin = reduce(hcat, remade_solution.u)
    sol_fin = sum(sol_fin, dims=1)


    if_display(
        Plots.scatter(
            data[1, :],
            data[2, :],
            xlabel="Time",
            ylabel="Arb. Units",
            label=["Data " nothing],
            markersize=2,
            color=:black,
            title=string(label_exp, " ", name_well),
            guidefontsize=guidefontsize,
            tickfontsize=tickfontsize,
            legendfontsize=legendfontsize,
            size=(y_size,x_size),
        ),
    )
    if_display(
        Plots.plot!(
            remade_solution.t,
            sol_fin[1, 1:end],
            xlabel="Time",
            ylabel="Arb. Units",
            label=[string("Fitting ", model) nothing],
            c=:red,
            guidefontsize=guidefontsize,
            tickfontsize=tickfontsize,
            legendfontsize=legendfontsize,
            size=(y_size,x_size),
        ),
    )
    if save_plot
        png(string(path_to_plot, label_exp, "_", model, "_", name_well, ".png"))
    end

    
end

function plot_NL_fit(save_plots,if_display,data,results_specific_well)

    name_well =   results_specific_well[2]
    label_exp =   results_specific_well[1]
    
end

function plot_seg_fit(save_plots,if_display,   data,results_specific_well)

    name_well =   results_specific_well[2]
    label_exp =   results_specific_well[1]
    
end


function plot_seg_fit(save_plots,if_display,   data,results_specific_well)

    name_well =   results_specific_well[2]
    label_exp =   results_specific_well[1]
    
end



path_to_data = "/Users/fabrizio.angaroni/Documents/JMAKi_utilities/real_dataset_tests/dataset/Monod_AA_detection/exp_s7/channel_1.csv"
path_to_annotation ="/Users/fabrizio.angaroni/Documents/JMAKi_utilities/real_dataset_tests/dataset/Monod_AA_detection/exp_s7/annotation.csv"



plot_data(
    "testing_plots", #label of the experiment
    path_to_data; # path to the folder to analyze
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
    display_plots=true,# display plots in julia or not
    save_plots=false, # save the plot or not
    overlay_plots=true, # true a single plot for all dataset false one plot per well
    avg_replicate=false, # if true the average between replicates
    x_size=300,
    do_blank_subtraction="avg_subtraction", # string on how to use blank (NO,avg_subtraction,time_avg)
    y_size =500,
    guidefontsize=18,
    tickfontsize=16,
    legendfontsize=10
    )

fit_log_lin = fit_one_file_Log_Lin(
    " ", #label of the experiment
    path_to_data; # path to the folder to analyze
    path_to_annotation=path_to_annotation,# path to the annotation of the wells
    avg_replicate=false, # if true the average between replicates is fitted. If false all replicate are fitted indipendelitly
    do_blank_subtraction="avg_subtraction", # string on how to use blank (NO,avg_subtraction,time_avg)

    )


plot_fit_of_file(
    fit_log_lin,
    path_to_data; # path to the folder to analyze
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
    path_to_plot="NA", # path where to save Plots
    do_blank_subtraction="avg_subtraction", # string on how to use blank (NO,avg_subtraction,time_avg)
    display_plots=true,# display plots in julia or not
    save_plots=false, # save the plot or not
    x_size=400,
    y_size =700,
    guidefontsize=18,
    tickfontsize=16,
    legendfontsize=10,
    pt_smoothing_derivative = 7 ,
)
  



model = "baranyi_richards"
 
lb_param = [0.001,0.1,0.0,0.01]
 ub_param =[0.1,5.0 ,1000.0,5.01]
param_guess =[0.01,1.0 ,500.0,1.01]
    
fit_od = fit_file_ODE(
        "test", #label of the experiment
        path_to_data, # path to the folder to analyze
        model, # string of the used model
        param_guess;
        path_to_annotation=path_to_annotation,# path to the annotation of the wells
        integrator=Tsit5(), # selection of sciml integrator
        lb = lb_param,
        ub =ub_param
    )
    


    plot_fit_of_file(
        fit_log_lin,
        path_to_data; # path to the folder to analyze
        path_to_annotation = path_to_annotation,# path to the annotation of the wells
        path_to_plot="NA", # path where to save Plots
        do_blank_subtraction="avg_subtraction", # string on how to use blank (NO,avg_subtraction,time_avg)
        display_plots=true,# display plots in julia or not
        save_plots=false, # save the plot or not
        x_size=400,
        y_size =700,
        guidefontsize=18,
        tickfontsize=16,
        legendfontsize=10,
        pt_smoothing_derivative = 7 ,
    )
      