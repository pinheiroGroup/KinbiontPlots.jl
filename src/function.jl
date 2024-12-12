
using CSV
using Plots
using Kinbiont

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
    correct_negative="remove", # if "thr_correction" it put a thr on the minimum value of the data with blank subracted, if "blank_correction" uses blank distrib to impute negative values
    thr_negative=0.01 ,
    blank_value = 0.0,
    blank_array = [0.0],
    x_size=400,
    y_size =750,
    guidefontsize=18,
    tickfontsize=16,
    legendfontsize=10
    )
        # creating the folder of data if one wants to save
        if save_plots == true
            mkpath(path_to_plot)
        end
  
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


    for well_name in names_of_cols[2:end]
        well_name = string(well_name)

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
                        label=[well_name],
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
                        label=[well_name],
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
                    title=string(label_exp, " ", well_name),
                    guidefontsize=guidefontsize,
                    tickfontsize=tickfontsize,
                    legendfontsize=legendfontsize,
                    size=(y_size,x_size),

                ),
            )
            if save_plots
                png(string(path_to_plot, label_exp, "_", well_name, ".png"))
            end
        end

    end
end







"""

    plot_fit_of_file(
    Kinbiont_results::Any;
    path_to_plot="NA", # path where to save Plots
    display_plots=true,# display plots in julia or not
    save_plots=false, # save the plot or not
    x_size=500,
    pt_smoothing_derivative = 7,
    y_size =750,
    guidefontsize=18,
    tickfontsize=16,
    legendfontsize=10,
    )

This function plot all the data from .csv file, note that assume that the first colums is the time

# Arguments:
- `Kinbiont_results`: The data struct of the results of a "fit  one .csv function" of Kinbiont


# Key Arguments:
-  `path_to_annotation::Any = missing`: The path to the .csv of annotation 
- `path_to_plot= "NA"`:String, path to save the plots.
-  `save_plot=false` :Bool, save the plot or not
- ` display_plots=true`:Bool,  Whether or not diplay the plot in julia
- `average_replicate=false` Bool, perform or not the average of replicates. Works only if an annotation path is provided
- `x_size=300` : X Size of the plot,
- `y_size=300` : y Size of the plot,
- `tickfontsize = 10` : size of the text of ticks in the plots
- `guidefontsize = 10` : size of the text of guide in the plots
- `legendfontsize = 10` : size of the text of legend in the plots

# Output:
- For this function the output are saved or displayed depending on the values of key arguments.
"""
function plot_fit_of_file(
       Kinbiont_results::Any;
        path_to_plot="NA", # path where to save Plots
        display_plots=true,# display plots in julia or not
        save_plots=false, # save the plot or not
        x_size=500,
        pt_smoothing_derivative = 7,
        y_size =750,
        guidefontsize=18,
        tickfontsize=16,
        legendfontsize=10,
        )
          # creating the folder of data if one wants to save
    if save_plots == true
        mkpath(path_to_plot)
    end
    # plotting standard fits
    if length(Kinbiont_results) == 4
       Kinbiont_results_matrix =Kinbiont_results[2]
       Kinbiont_fits =Kinbiont_results[3]
       Kinbiont_data =Kinbiont_results[4]



        if display_plots
            if_display = display
        else
            if_display = identity
        end

      for i in eachindex(Kinbiont_fits)
          fit_temp =Kinbiont_fits[i]
          data_temp =Kinbiont_data[i]
          y_data_temp =data_temp[2,:]
          x_data_temp =data_temp[1,:]
          model_string =Kinbiont_results_matrix[3,i+1]
          well_name =Kinbiont_results_matrix[2,i+1]
          label_exp =Kinbiont_results_matrix[1,i+1]

          if_display(
            Plots.scatter(
                x_data_temp,
                y_data_temp,
                xlabel="Time",
                ylabel="Arb. Units",
                label=["Data " nothing],
                markersize=2,
                color=:black,
                title=string(label_exp, " ", well_name),
                guidefontsize=guidefontsize,
                tickfontsize=tickfontsize,
                legendfontsize=legendfontsize,
                size=(y_size,x_size),
            ),
            )
            if_display(
                Plots.plot!(
                    x_data_temp,
                    fit_temp,
                    xlabel="Time",
                    ylabel="Arb. Units",
                    label=[string("Fitting ", model_string) nothing],
                    c=:red,
                    guidefontsize=guidefontsize,
                    tickfontsize=tickfontsize,
                    legendfontsize=legendfontsize,
                    size=(y_size,x_size),
                ),
            )
            if save_plots
                png(string(path_to_plot, label_exp, "_", model_string, "_", well_name, ".png"))
            end



        end
   

        
    elseif Kinbiont_results[1] == "Log-Lin" 

       Kinbiont_results_matrix =Kinbiont_results[2]
       Kinbiont_fits =Kinbiont_results[3]
       Kinbiont_data =Kinbiont_results[4]
       Kinbiont_confidence_intervals =Kinbiont_results[5]

        if display_plots
            if_display = display
        else
            if_display = identity
        end

      for i in eachindex(Kinbiont_fits)


        if !ismissing(Kinbiont_fits[i][1])
          fit_temp =Kinbiont_fits[i]
          data_temp =Kinbiont_data[i]
          temp_ci =Kinbiont_confidence_intervals[i]
          y_fit_temp =fit_temp[:,2]
          x_fit_temp =fit_temp[:,1]
          y_data_temp =data_temp[2,:]
          x_data_temp =data_temp[1,:]
          well_name =Kinbiont_results_matrix[2,i+1]
          label_exp =Kinbiont_results_matrix[1,i+1]
          coeff_1 =Kinbiont_results_matrix[12,i+1]
          coeff_2 =Kinbiont_results_matrix[7,i+1]


        y_fit_temp =fit_temp[:,2]
        x_fit_temp =fit_temp[:,1]
        y_data_temp =data_temp[2,:]
        x_data_temp =data_temp[1,:]


        if_display(
            Plots.scatter(
                x_data_temp,
                log.(y_data_temp),
                xlabel="Time",
                ylabel="Log(Arb. Units)",
                label=["Log Data " nothing],
                markersize=1,
                color=:black,
                title=string(label_exp, " ", well_name),
                guidefontsize=guidefontsize,
                tickfontsize=tickfontsize,
                legendfontsize=legendfontsize,
                size=(y_size,x_size),
    
            ),
        )
    
        if_display(
            Plots.plot!(
                x_fit_temp,
                y_fit_temp,
                ribbon=temp_ci,
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
                [x_fit_temp[1], x_fit_temp[end]],
                c=:black,
                label=[string("Window of exp. phase ") nothing],
                guidefontsize=guidefontsize,
                tickfontsize=tickfontsize,
                legendfontsize=legendfontsize,
                size=(y_size,x_size),
    
            ),
        )
        if save_plots
            png(string(path_to_plot, label_exp, "_Log_Lin_Fit_", well_name, ".png"))
        end
  
        N0 = exp.(coeff_1)
        Theoretical_fitting_exp = N0.* exp.(x_fit_temp .* coeff_2)
    
        if_display(
            Plots.scatter(
                x_data_temp,
                y_data_temp,
                xlabel="Time",
                ylabel="Arb. Units",
                label=["Data " nothing],
                markersize=1,
                color=:black,
                title=string(label_exp, " ", well_name),
                guidefontsize=guidefontsize,
                tickfontsize=tickfontsize,
                legendfontsize=legendfontsize,
                size=(y_size,x_size),
    
            ),
        )
    
        if_display(
            Plots.plot!(
                x_fit_temp,
                Theoretical_fitting_exp,
                xlabel="Time ",
                ylabel="Arb. Units",
                label=[string("Fitting Exp. ") nothing],
                c=:red,
                guidefontsize=guidefontsize,
                tickfontsize=tickfontsize,
                legendfontsize=legendfontsize,
                size=(y_size,x_size),
    
            ),
        )
        if_display(
            Plots.vline!(
                [x_fit_temp[1], x_fit_temp[end]],
                c=:black,
                label=[string("Window of exp. phase ") nothing],
                guidefontsize=guidefontsize,
                tickfontsize=tickfontsize,
                legendfontsize=legendfontsize,
                size=(y_size,x_size),
    
            ),
        )
        if save_plots
            png(string(path_to_plot, label_exp, "_exp_Fit_", well_name, ".png"))
        end
        data = Matrix(transpose(hcat(x_data_temp,y_data_temp)))



        specific_gr = Kinbiont.specific_gr_evaluation(data, pt_smoothing_derivative)
    
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
                title=string(label_exp, " ", well_name),
                guidefontsize=guidefontsize,
                tickfontsize=tickfontsize,
                legendfontsize=legendfontsize,
                size=(y_size,x_size),
    
            ),
        )
        if_display(
            Plots.vline!(
                [x_fit_temp[1], x_fit_temp[end]],
                c=:black,
                label=[string("Window of exp. phase ") nothing],
                   guidefontsize=guidefontsize,
            tickfontsize=tickfontsize,
            legendfontsize=legendfontsize,
            size=(y_size,x_size),
        ) )
        if save_plots
            png(string(path_to_plot, label_exp, "_specific_gr_dynamics_", well_name, ".png"))
        end

    end
end

    else    # plotting segmentation fits

       Kinbiont_results_matrix =Kinbiont_results[2]
       Kinbiont_data =Kinbiont_results[4]
       Kinbiont_cp_intervals =Kinbiont_results[3]
        well_names = unique(Kinbiont_results_matrix[2,2:end])


        if display_plots
            if_display = display
        else
            if_display = identity
        end

      for i in eachindex(Kinbiont_fits)
          fit_temp =Kinbiont_fits[i]
          data_temp =Kinbiont_data[i]
          temp_cp =Kinbiont_cp_intervals[i]
          y_fit_temp =fit_temp[:,2]
          x_fit_temp =fit_temp[:,1]
          y_data_temp =data_temp[2,:]
          x_data_temp =data_temp[1,:]
          model_string =Kinbiont_results_matrix[3,i+1]
          well_name = well_names[i]
          label_exp =Kinbiont_results_matrix[1,i+1]

          if_display(
            Plots.scatter(
                x_data_temp,
                y_data_temp,
                xlabel="Time",
                ylabel="Arb. Units",
                label=["Data " nothing],
                markersize=2,
                color=:black,
                title=string(label_exp, " ", well_name),
                guidefontsize=guidefontsize,
                tickfontsize=tickfontsize,
                legendfontsize=legendfontsize,
                size=(y_size,x_size),
            ),
            )



            if_display( Plots.vline!(
                [temp_cp],
                c=:black,
                label=[string("Change point") nothing],
                guidefontsize=guidefontsize,
                tickfontsize=tickfontsize,
                legendfontsize=legendfontsize,
                size=(y_size,x_size),
                ) )
            if save_plots
                png(string(path_to_plot, label_exp, "_segmented_fit_", well_name, ".png"))
            end



        end
   

        
        

    end


    
      
end    


export plot_data
export plot_fit_of_file