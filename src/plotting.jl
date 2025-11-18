using Plots 
using TOML 
using Dates 

include(joinpath(@__DIR__, "constants.jl")) 

#input_file_path = joinpath(@__DIR__, "..", "main", "input_values.toml") 
#input_dict = load_parameters(input_file_path) 


function PlotTrajectory( sol, input ) 

    HeliopauseLine = input["distance_HP"]
    params = input["p"]
    plot_color = input["plot_color"]
    plane = input["plane"] 
    plot_sun = input["plot_sun"]
    plot_vectors = input["plot_vectors"]
    n_grid = input["n_grid"]
    plot_color = input["plot_color"]
    dist_measure = input["dist_measure"] 
    plasma_model = input["plasma_model"] 

    # --- Extract components ---
    sol_time = sol.t ./ yr
    sol_x = [u[1] for u in sol.u] ./ AU
    sol_y = [u[2] for u in sol.u] ./ AU
    sol_z = [u[3] for u in sol.u] ./ AU
    sol_speed = [norm(u[4:6]) for u in sol.u] ./ 1e3   # in km/s

    # --- Downsample trajectory if too long ---
    npts = length(sol_time)
    if npts > 1000
        idx = round.(Int, range(1, npts, length=1000))
        sol_time = sol_time[idx]
        sol_x, sol_y, sol_z, sol_speed = sol_x[idx], sol_y[idx], sol_z[idx], sol_speed[idx]
    end

    # --- Helper function to make a single plot ---
    function make_plane_plot(xdata, ydata, xlabel_, ylabel_, title_plane; is_xz=false)
        xleft, xright = (maximum([70, minimum(xdata)-9]), minimum([130, maximum(xdata)+9])) #(70, 130)
        ybottom, ytop = (minimum(ydata)-12, maximum(ydata)+12)
        boundary_scale = 40.0

        # Determine color mapping (time vs speed)
        if plot_color == "time"
            sol_coloring = sol_time
            cbar_title = "Time [yr]"
        elseif plot_color == "speed"
            sol_coloring = sol_speed
            cbar_title = "Speed [km/s]"
        else
            sol_coloring = sol_time
            cbar_title = "Time [yr]"
        end

        trajectory_duration = (input["max_time"] - input["min_time"]) / yr 
        title_text =
            "Beta: $(params[1]); Qm: $(params[2]); Plasma Model: $(plasma_model); \n" *
            "Duration: $(trajectory_duration) yrs; Plane: $(title_plane)" 


        plt = plot(
            xdata, ydata,
            linez = sol_coloring,
            linewidth = 1.5,
            color = :plasma,
            colorbar = :true,
            colorbar_title = cbar_title,
            xlabel = xlabel_,
            ylabel = ylabel_,
            title = title_text, #"Trajectory: $title_plane",
            size = (900, 550),
            legend = false,
            grid = true,
            gridlinewidth = 2,
            margin = 2Plots.mm,
            alpha = 0.95,
            xlims = (xleft, xright), #(maximum([xleft, minimum(xdata)-5]), minimum([xright, maximum(xdata)+5])), 
            ylims = (ybottom, ytop) #(minimum(ydata)-12, maximum(ydata)+12)
        )

        # --- Plot plasma velocity arrows ---
        if plot_vectors
            xgrid = range(xleft, xright, length=n_grid)
            #ygrid = range(minimum([minimum(ydata),-boundary_scale]), 
            #                maximum([maximum(ydata),boundary_scale]), length=n_grid) 
            ygrid = range(ybottom, ytop, length=n_grid) 

            X, Y, U, V = Float64[], Float64[], Float64[], Float64[]
            for x in xgrid, y in ygrid
                uvec = is_xz ? [x*AU, 0.0, y*AU, 0, 0, 0] : [x*AU, y*AU, 0.0, 0, 0, 0]
                v = PlasmaVelocity(uvec, input)
                push!(X, x); push!(Y, y)
                if is_xz
                    push!(U, v[1]/1e4 / AU)
                    push!(V, v[3]/1e4 / AU)
                else
                    push!(U, v[1]/1e4 / AU)
                    push!(V, v[2]/1e4 / AU)
                end
            end

            quiver!(plt, X, Y, quiver=(U, V), color=:gray, alpha=0.8, linewidth=0.8, label="Plasma flow")
        end

        # --- Plot Sun ---
        if plot_sun
            scatter!(plt, [0.0], [0.0], markersize=8, color=:yellow, label="Sun")
        end

        # --- Draw heliopause boundary ---
        if HeliopauseLine != 0.0
            if dist_measure == "flat"
                plot!(plt, [HeliopauseLine/AU, HeliopauseLine/AU], [-boundary_scale, boundary_scale],
                      color=:red, linewidth=2, linestyle=:solid, label="Heliopause")
            elseif dist_measure == "spherical"
                y_arc = range(-boundary_scale, boundary_scale, length=200)
                x_arc = sqrt.((HeliopauseLine / AU).^2 .- y_arc.^2)
                plot!(plt, x_arc, y_arc, color=:red, linewidth=2, linestyle=:solid, label="Heliopause")
            end
        end

        # --- Magnetic Field Visualization ---
        y_pos = ydata[1] + 5.0 #25.0
        u_inside  = [HeliopauseLine - 1.0 * AU, 0.0, 0.0, 0, 0, 0]
        u_outside = [HeliopauseLine + 1.0 * AU, 0.0, 0.0, 0, 0, 0]

        B_in, B_out = B_field(u_inside, input), B_field(u_outside, input)
        Bmag_in, Bmag_out = norm(B_in), norm(B_out)

        annotate!(plt, (HeliopauseLine/AU - 7.0, y_pos - 4.0, text("$(round(Bmag_in*1e9, digits=2)) nT", 8, :blue)))
        annotate!(plt, (HeliopauseLine/AU + 7.0, y_pos - 4.0, text("$(round(Bmag_out*1e9, digits=2)) nT", 8, :green)))
        annotate!(plt, (HeliopauseLine/AU - 7.0, y_pos + 6.0, text("Heliosphere", 8, :blue)))
        annotate!(plt, (HeliopauseLine/AU + 7.0, y_pos + 6.0, text("ISM", 8, :green)))        

        if is_xz
            # --- XZ plane → symbol indicators ---
            inside_symbol = B_in[2] > 0 ? :x : :circle
            outside_symbol = B_out[2] > 0 ? :x : :circle

            # Inside
            if inside_symbol == :x
                scatter!(plt, [HeliopauseLine/AU - 7.0], [y_pos],
                        marker=(inside_symbol, 8),
                        color=:blue,
                        #markerstrokewidth=1.5,
                        label="")
            else
                scatter!(plt, [HeliopauseLine/AU - 7.0], [y_pos],
                        marker=(inside_symbol, 8),
                        markercolor=:white,
                        markerstrokecolor=:blue,
                        markerstrokewidth=1.5,
                        label="")
            end

            # Outside
            if outside_symbol == :x
                scatter!(plt, [HeliopauseLine/AU + 7.0], [y_pos],
                        marker=(outside_symbol, 8),
                        color=:green,
                        #markerstrokewidth=1.5,
                        label="")
            else
                scatter!(plt, [HeliopauseLine/AU + 7.0], [y_pos],
                        marker=(outside_symbol, 8),
                        markercolor=:white,
                        markerstrokecolor=:green,
                        markerstrokewidth=1.5,
                        label="")
            end
        #end


        else
            # --- XY plane → draw field vectors ---
            x_B_in, x_B_out = [[HeliopauseLine/AU - 5.0], [HeliopauseLine/AU + 5.0]]
            y_B_in, y_B_out = [[y_pos], [y_pos]]
            U_B_in, U_B_out = [[B_in[1]/Bmag_in].* 5, [B_out[1]/Bmag_out].* 5]  
            V_B_in, V_B_out = [[B_in[2]/Bmag_in].* 5, [B_out[2]/Bmag_out].* 5] 

            quiver!(plt, x_B_in, y_B_in, quiver=(U_B_in, V_B_in),
                    color=:blue, alpha=0.8, linewidth=1.5, label="Magnetic field (inside)")            
            quiver!(plt, x_B_out, y_B_out, quiver=(U_B_out, V_B_out),
                    color=:green, alpha=0.8, linewidth=1.5, label="Magnetic field (outside)") 
        end

        # --- Final position marker and Q/m label ---
        scatter!(plt, [xdata[end]], [ydata[end]], color=:black, markersize=5, label="Final position")
        if params !== nothing
            annotate!(plt, (xdata[end]-2, ydata[end]-4, text("Q/m = $(round(params[2], sigdigits=4))", 8, :black)))
        end

        return plt
    end

    # --- Select projection plane ---
    if plane == "xz"
        plt = make_plane_plot(sol_x, sol_z, "x [AU]", "z [AU]", "x–z plane"; is_xz=true)
    elseif plane == "xy"
        plt = make_plane_plot(sol_x, sol_y, "x [AU]", "y [AU]", "x–y plane"; is_xz=false)
    elseif plane == "both"
        plt_xy = make_plane_plot(sol_x, sol_y, "x [AU]", "y [AU]", "x–y plane"; is_xz=false)
        plt_xz = make_plane_plot(sol_x, sol_z, "x [AU]", "z [AU]", "x–z plane"; is_xz=true)
        plt = plot(plt_xy, plt_xz, layout=(2,1), size=(950,900), margin=4Plots.mm)
    else
        error("Invalid plane argument. Use \"xz\", \"xy\", or \"both\".")
    end

    # Save figure if requested 
    if get(input, "save_fig", false) == true
        q_over_m = params[2]   # Q/m stored in p = [something, q/m, ...]

        # Timestamp: YYYY-MM-DD-HH-MM
        date_str = Dates.format(Dates.now(), "yyyy-mm-dd--HH-MM-SS")

        # Build filename: Year-Month-Day-Hour-Min-q_over_m.png
        fname = "$(date_str)--Qm_$(q_over_m)--$(plasma_model).png"

        # Save under /data relative to project root
        filename = joinpath(@__DIR__, "..", "data", fname)

        try
            savefig(plt, filename)
            println("Saved figure to: $filename")
        catch e
            @warn "Failed to save figure" error=e
        end
    end


    return plt 
end 
