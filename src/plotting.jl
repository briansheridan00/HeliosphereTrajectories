using Plots 
using TOML 
using Dates 

include(joinpath(@__DIR__, "constants.jl"))  


# Function to save the plot 
function save_plot_if_requested(plt, input, params)
    if get(input, "save_fig", false)
        q_over_m = params[2]
        date_str = Dates.format(Dates.now(), "yyyy-mm-dd--HH-MM-SS")
        fname = "$(date_str)--Qm_$(q_over_m)--$(input["plasma_model"])--$(input["B_model"]).png"
        filename = joinpath(@__DIR__, "..", "data", fname)

        try
            savefig(plt, filename)
            println("Saved figure to: $filename")
        catch e
            @warn "Failed to save figure" error=e
        end
    end
end


### --- Main plotting function --- ###
function Plotter(sol, input; plot_B_fields=true, charges=false) 
    params = input["p"] 

    # Base trajectory plot
    plt_base = SlicePlot(sol, input)

    # Case 1: No magnetic field plots requested 
    if !plot_B_fields
        if charges != false
            qm_plot = qm_over_time(sol, charges, input["charging_type"])
            layout_def = @layout [a; b]
            final_plot = plot(
                plt_base,
                qm_plot;
                layout=layout_def,
                size=(1200, 650),
                margin=4Plots.mm
            )
            save_plot_if_requested(final_plot, input, params)
            return final_plot
        else
            save_plot_if_requested(plt_base, input, params)
            return plt_base
        end
    end

    # Case 2: Magnetic field plots ARE requested

    # Time values
    t_start = input["min_time"]
    t_end   = input["max_time"]
    t_mid   = (t_end - t_start) / 2.0
    t_other = 4.5 * (t_end - t_start) / 5.0
    nx, nt, nza = 150, 170, 25 #180, 200, 30

    # Magnetic field plots
    plots = [
        plt_base,
        plot_B_polarity(input; radius=90, time=t_start, n_x=nx, n_t=nt, n_z=nza, visual_mode="xz", msize=4.5),
        plot_B_polarity(input; radius=90, time=t_mid,   n_x=nx, n_t=nt, n_z=nza, visual_mode="xz", msize=4.5),
        plot_B_polarity(input; radius=50, time=t_other, n_x=nx, n_t=nt, n_z=nza, visual_mode="xz", msize=4.5)
    ]

    # Add Q/m plot if requested
    if charges != false
        qm_plot = qm_over_time(sol, charges, input["charging_type"])

        layout_def = @layout [
            a{0.7w} [grid(3,1)]
            b{0.3h}
        ]

        final_plot = plot(
            plots...,
            qm_plot;
            layout=layout_def,
            size=(1800, 1200),
            margin=7Plots.mm
        )

        save_plot_if_requested(final_plot, input, params)
        return final_plot
    end

    # If no Q/m plot is added, return only B-field composite
    layout_def = @layout [
        a{0.7w} [grid(3,1)]
    ]

    final_plot = plot(
        plots...;
        layout=layout_def,
        size=(1400, 600),
        margin=4Plots.mm
    )

    save_plot_if_requested(final_plot, input, params)
    return final_plot
end



function PlotterOld(sol, input; plot_B_fields=true, charges=false) 
    params = input["p"] 

    # First base plot
    plt_base = SlicePlot(sol, input)
 
    if plot_B_fields
            # Decide what "end" means → final time of solution
            t_start = input["min_time"]
            t_end = input["max_time"]
            t_mid = (t_end - t_start) / 2.0 
            t_other = 4.5 * (t_end - t_start) / 5.0 
            t_other2 = 1.0 * (t_end - t_start) / 22.0 
            nx, nt, nza, nzb = [180, 200, 30, 30]

            plots = [
                plt_base,
                plot_B_polarity(input; radius=90, time=t_start, n_x=nx, n_t=nt, n_z=nza, visual_mode="xz", msize=4.5),
                plot_B_polarity(input; radius=90, time=t_mid, n_x=nx, n_t=nt, n_z=nza, visual_mode="xz", msize=4.5),
                plot_B_polarity(input; radius=50, time=t_other , n_x=nx, n_t=nt, n_z=nza, visual_mode="xz", msize=4.5)
            ] 

            layout_def = @layout [
                a{0.7w} [grid(3,1)]
            ]

            final_plot = plot(
                plots...;
                layout=layout_def,
                size=(1400, 600),
                margin=4Plots.mm
            )

        # Save plot if asked
        save_plot_if_requested(final_plot, input, params)
        return final_plot
    end
    save_plot_if_requested(plt_base, input, params)
    return plt_base
end



# --- Create a plot of the plane of interest --- #
function SlicePlot(sol, input)  
    plane = "xz" #input["plane"] 

    # --- Extract and downsample trajectory ---
    sol_time = sol.t ./ yr
    sol_x = [u[1] for u in sol.u] ./ AU
    sol_y = [u[2] for u in sol.u] ./ AU
    sol_z = [u[3] for u in sol.u] ./ AU
    trajectory_duration = (sol_time[end] - sol_time[1])  # in years 
    sol_speed = [norm(u[4:6]) for u in sol.u] ./ 1e3   # in km/s   
    npts = length(sol_time)
    if npts > 1000
        idx = round.(Int, range(1, npts, length=1000))
        sol_time = sol_time[idx]
        sol_x, sol_y, sol_z, sol_speed = sol_x[idx], sol_y[idx], sol_z[idx], sol_speed[idx]
    end     

    # --- Define the abcissa and ordinate axes data and limits --- 
    xdata = sol_x 
    ydata = plane == "xz" ? sol_z : sol_y 
    x_label = "x [AU]"
    y_label = plane == "xz" ? "z [AU]" : "y [AU]" 
    title_text = "Trajectory ---> Beta: $(input["beta_value"]); Charging: $(input["charging_type"]);\n" * 
                "Plasma: $(input["plasma_model"]), B_field: $(input["B_model"]); Duration: $(round(trajectory_duration,digits=3)) yr" 

    input["xleft"] =  minimum( [70, minimum(xdata)-9] )
    input["xright"] = minimum( [130, maximum(xdata)+9] ) 
    input["ybottom"] =  minimum(ydata)-12
    input["ytop"] = maximum(ydata)+12                

    xleft, xright = ( minimum( [70, minimum(xdata)-9] ), minimum( [130, maximum(xdata)+9] ) ) #(70, 120)
    ybottom, ytop = ( minimum(ydata)-12, maximum(ydata)+12 ) #(-90, 90) 

    # --- Define and return the plot --- 
    plt = Plots.plot(
        xdata, ydata, linez = sol_speed, linewidth = 1.5, color = :plasma, colorbar = :true, colorbar_title = "Speed [km/s]",
        xlabel = x_label, ylabel = y_label, title = title_text, size = (900, 600), legend = false, grid = true,
        gridlinewidth = 2, margin = 4Plots.mm, alpha = 0.95,
        xlims = (xleft, xright),  ylims = (ybottom, ytop), 
        titlefont=10, guidefont=8, tickfont=8 )

    # --- Final position marker and Q/m label ---
    scatter!(plt, [xdata[end]], [ydata[end]], color=:black, markersize=5) 

    # --- Plot Plasma Field --- 
    is_xz =  plane == "xz" ? true : false 
    Plasma!(plt, input; plot_vectors=input["plot_vectors"], n_grid=20, is_xz=is_xz)  

    # --- Plot the boundary lines --- 
    BoundaryLines!(plt, input; boundary_scale=0.95)

    # --- Plot the sun --- 
    scatter!(plt, [0], [0], color=:yellow, markerstrokecolor=:black, markersize=5)

    # --- Annotate time steps along trajectory --- 
    AnnotateTimes!(plt, input, sol_time, xdata, ydata; n_labels=7)

    return plt  
end 


# --- Annotate the times along the trajectory --- 
function AnnotateTimes!(plt, input, sol_time, xdata, ydata; n_labels=7)
    if input["annotate_times"] == true 
        idxs = round.(Int, range(1, length(sol_time), length=n_labels))
        for i in idxs
            annotate!(plt, (xdata[i], ydata[i]-1.0, text("t: $(round(sol_time[i], digits=2)) yr", 7, :black))) 
        end 
    end 
    return plt 
end 


# --- Plot the plasma velocity field --- # 
function Plasma!(plt, input; plot_vectors=true, n_grid=20, is_xz=false)
    if !plot_vectors
        return plt
    end
    # Bounds 
    xgrid = range(input["xleft"], input["xright"], length=n_grid)
    ygrid = range(input["ybottom"], input["ytop"], length=n_grid)
    X = Float64[]; Y = Float64[]; U = Float64[]; V = Float64[]
    for x in xgrid, y in ygrid
        uvec = is_xz ?
            [x*AU, 0.0, y*AU, 0.0, 0.0, 0.0] :
            [x*AU, y*AU, 0.0, 0.0, 0.0, 0.0]
        v = PlasmaVelocity(uvec, input)
        push!(X, x)
        push!(Y, y)
        if is_xz
            push!(U, v[1]/1e4 / AU)
            push!(V, v[3]/1e4 / AU)
        else
            push!(U, v[1]/1e4 / AU)
            push!(V, v[2]/1e4 / AU)
        end
    end
    quiver!( plt, X, Y, quiver=(U, V), color=:gray, alpha=0.8, linewidth=0.8, label="Plasma flow" ) 
    return plt
end


# --- Plot the boundary lines --- # 
function BoundaryLines!(plt, input; boundary_scale=0.95)
    dist_measure = input["dist_measure"]
    HeliopauseLine = input["distance_HP"] /AU
    TerminationLine = input["distance_TS"] /AU

    # --- Draw Heliopause Boundary ---
    if HeliopauseLine != 0.0
        ylim_HP = boundary_scale * HeliopauseLine
        if dist_measure == "flat"
            plot!(plt, [HeliopauseLine, HeliopauseLine], [-ylim_HP, ylim_HP],
                    color=:red, linewidth=2, linestyle=:solid, label="Heliopause")
        elseif dist_measure == "spherical"
            y_arc = range(-ylim_HP, ylim_HP, length=200)
            x_arc = sqrt.((HeliopauseLine).^2 .- y_arc.^2)
            plot!(plt, x_arc, y_arc, color=:red, linewidth=2, linestyle=:solid, label="Heliopause")
        end
    end

    # --- Draw Termination Shock Line --- 
    if TerminationLine != 0.0
        ylim_TS = boundary_scale * TerminationLine
        if dist_measure == "flat"
            plot!(plt, [TerminationLine, TerminationLine], [-ylim_TS, ylim_TS],
                    color=:green, linewidth=2, linestyle=:solid, label="Termination")
        elseif dist_measure == "spherical"
            y_arc = range(-ylim_TS, ylim_TS, length=200)
            x_arc = sqrt.((TerminationLine).^2 .- y_arc.^2)
            plot!(plt, x_arc, y_arc, color=:green, linewidth=2, linestyle=:solid, label="Termination")
        end
    end 

    return plt 
end     



### --- Plot the magnetic fields over time and space --- ### 

# --- Helper functions --- 
flat_polarity(B_vector) = sign(B_vector[2])

spherical_polarity(u, B_vector) = begin
    phi = atan(u[2], u[1])
    e_phi = [-sin(phi), cos(phi), 0.0]
    sign(dot(B_vector, e_phi))
end

# --- Calculate magnetic field polarity over regimes of interest --- 
function B_data(input; radius=90.0, time=0.0, n_x=20, n_t=20, n_z=20, visual_mode="tz")

    @assert haskey(input, "min_time") && haskey(input, "max_time")

    x_span_low, x_span_high = (35, 110)
    y_span_value = 35

    if visual_mode == "tz"
        xspan = range(input["min_time"], input["max_time"], n_t) #time axis
        yspan = range(-y_span_value, y_span_value, n_z) #z axis
    elseif visual_mode == "xz"
        xspan = range(x_span_low, x_span_high, n_x) #x axis
        #println("xspan[end]: ", xspan[end])
        yspan = range(-y_span_value, y_span_value, n_z) #z axis 
    else 
        error("visual_mode not recognised")
    end 

    B_pol = zeros(length(xspan), length(yspan))

    if visual_mode == "tz"
        for (i, x) in enumerate(xspan)
            for (j, y) in enumerate(yspan)
                #println("\n x=$(x) and y=$(y)")
                u = [radius*AU, 0.0, y*AU, 0.0, 0.0, 0.0]
                B_vector = B_field(u, input; t=x)
                #println("B_vector: $(B_vector)")
                if input["dist_measure"] == "flat"
                    pol = flat_polarity(B_vector)
                    B_pol[i,j] = pol
                    #println("B pol: $(pol)")
                else 
                    pol = spherical_polarity(u, B_vector)
                    B_pol[i,j] = pol 
                    #println("B pol: $(pol)")
                end 
            end
        end
    elseif visual_mode == "xz"
        for (i, x) in enumerate(xspan)
            for (j, y) in enumerate(yspan)
                #println("\n x=$(x) and y=$(y)")
                u = [x*AU, 0.0, y*AU, 0.0, 0.0, 0.0]
                B_vector = B_field(u, input; t=time)
                #println("B_vector: $(B_vector)")
                if input["dist_measure"] == "flat"
                    pol = flat_polarity(B_vector)
                    B_pol[i,j] = pol
                    #println("B pol: $(pol)")
                else 
                    pol = spherical_polarity(u, B_vector)
                    B_pol[i,j] = pol 
                    #println("B pol: $(pol)")
                end 
            end
        end
    end 
    return xspan, yspan, B_pol 
end

# --- Plot the magnetic field plot --- # 
function plot_B_polarity(input; radius=90.0, time="end", n_x=20, n_t=20, n_z=20, visual_mode="tz", msize=2.0)
    if time == "end"
        time_val = Float64(input["max_time"])
    elseif time == "start" 
        time_val = Float64(input["min_time"]) 
    elseif typeof(time) <: Union{Float64, Float32, Int} 
        time_val = time 
    else 
        time_val = 0.0
    end 

    # Create the data # Outputs: time in seconds, x in AU, z in AU. 
    xvalues, yvalues, Bpolarities = B_data(input; radius=radius, time=time_val, n_x=n_x, n_t=n_t, n_z=n_z, visual_mode=visual_mode)

    # Flatten the grids 
    if visual_mode == "tz"
        xvec = repeat(xvalues ./ yr, outer=length(yvalues)) # outer loop 
        yvec = repeat(yvalues, inner=length(xvalues)) # inner loop 
        polvec = vec(Bpolarities)    
    elseif visual_mode == "xz"
        xvec = repeat(xvalues, outer=length(yvalues))
        yvec = repeat(yvalues, inner=length(xvalues))  
        polvec = vec(Bpolarities) 
    end 

    #println("Start and end of xvec: $(xvec[1]) and $(xvec[end])")

    # Labels 
    x_label = visual_mode == "tz" ? "Time [yr]" : "x [AU]"
    y_label = visual_mode == "tz" ? "z [AU]" : "z [AU]"  
    title_string = "B_field polarity - Model: $(input["B_model"]); "
    #println("vis mode: ", visual_mode)
    title_append = visual_mode == "tz" ? "Radius: $(round(radius, digits=2)) AU" : "Time $(round(time_val / yr, digits=2)) yr"
    title_string *= title_append 
    
    # Map polarity to colors: -1 → blue, +1 → red
    colors = map(ccc -> ccc == 1 ? :red : :blue, polvec) 
    
    plt = Plots.scatter(
        xvec, yvec, c = colors, ms=msize, #2.0, #3.5, 
        markerstrokecolor=:black, markerstrokewidth=0.5,
        xlabel=x_label, ylabel=y_label, title=title_string, legend=false, colorbar=:true, margin = 2Plots.mm, 
        titlefont=8, guidefont=6, tickfont=6
    ) 
    return plt 
end


# --- Charge to mass ratio over time --- 
function qm_over_time(trajectory_solution, trajectory_qm_values, charging_type)
    plt = plot(trajectory_solution.t ./ yr, trajectory_qm_values, 
            xlabel="Time (yr)", ylabel="Q/m",
            size=(800, 250), label="Q/m", margin=5Plots.mm, linewidth=2.5,
            title="Charge to mass ratio – $(charging_type) charging", 
            ylim=(-0.5, maximum(trajectory_qm_values)+1)) 

    return plt 
end 

