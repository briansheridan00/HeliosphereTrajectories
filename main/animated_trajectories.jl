using Pkg
using TOML 
using LinearAlgebra
using DifferentialEquations
using Random
using Plots
using Interpolations 

include(joinpath(@__DIR__, "..", "src", "utils.jl"))
include(joinpath(@__DIR__, "..", "src", "constants.jl"))
include(joinpath(@__DIR__, "..", "src", "magnetic_field.jl"))
include(joinpath(@__DIR__, "..", "src", "plasma_field.jl"))
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "charging.jl"))
include(joinpath(@__DIR__, "..", "src", "trajectory.jl")) 

using JLD2, FileIO
data = load(joinpath(@__DIR__, "..", "data", "charging_dict.jld2"))
charging_dict = data["charging_dict"]

input_file = joinpath(@__DIR__, "..", "main", "input_values.toml")
input_parameters = load_parameters(input_file) 

res_traj, saved_charges = ComputeTrajectory(input_parameters); 

sol_time = res_traj.t #./ yr
sol_x = [u[1] for u in res_traj.u] #./ AU
sol_y = [u[2] for u in res_traj.u] #./ AU
sol_z = [u[3] for u in res_traj.u] #./ AU
trajectory_duration =  (sol_time[end] - sol_time[1])   
sol_speed = [norm(u[4:6]) for u in res_traj.u] #./ 1e3   # in km/s   
npts = length(sol_time)
if npts > 1000
    idx = round.(Int, range(1, npts, length=1000))
    sol_time = sol_time[idx]
    sol_x, sol_y, sol_z, sol_speed = sol_x[idx], sol_y[idx], sol_z[idx], sol_speed[idx]
end 



# --- Helper functions --- 
flat_polarity(B_vector) = sign(B_vector[2])

spherical_polarity(u, B_vector) = begin
    phi = atan(u[2], u[1])
    e_phi = [-sin(phi), cos(phi), 0.0]
    sign(dot(B_vector, e_phi))
end

# --- Calculate magnetic field polarity over regimes of interest --- 
function B_data(input; radius=90.0, timestamp=0.0, n_x=20, n_t=20, n_z=20, visual_mode="tz", xspan_range=(70,110), yspan_range=20)

    @assert haskey(input, "min_time") && haskey(input, "max_time")

    x_span_low, x_span_high = xspan_range #(0.1, 110)
    y_span_value = yspan_range #55

    if visual_mode == "tz"
        xspan = range(timestamp + input["min_time"], timestamp + input["max_time"], n_t) #time axis
        yspan = range(-y_span_value, y_span_value, n_z) #z axis
    elseif visual_mode == "xz"
        xspan = range(x_span_low, x_span_high, n_x) #x axis 
        yspan = range(-y_span_value, y_span_value, n_z) #z axis 
    else 
        error("visual_mode not recognised")
    end 

    B_pol = zeros(length(xspan), length(yspan))

    if visual_mode == "tz"
        for (i, x) in enumerate(xspan)
            for (j, y) in enumerate(yspan) 
                u = [radius*AU, 0.0, y*AU, 0.0, 0.0, 0.0]
                B_vector = B_field(u, input; t=x) 
                if input["dist_measure"] == "flat"
                    pol = flat_polarity(B_vector)
                    B_pol[i,j] = pol 
                else 
                    pol = spherical_polarity(u, B_vector)
                    B_pol[i,j] = pol  
                end 
            end
        end
    elseif visual_mode == "xz"
        for (i, x) in enumerate(xspan)
            for (j, y) in enumerate(yspan) 
                u = [x*AU, 0.0, y*AU, 0.0, 0.0, 0.0]
                B_vector = B_field(u, input; t=timestamp)
                if input["dist_measure"] == "flat"
                    pol = flat_polarity(B_vector)
                    B_pol[i,j] = pol 
                else 
                    pol = spherical_polarity(u, B_vector)
                    B_pol[i,j] = pol  
                end 
            end
        end
    end 
    return xspan, yspan, B_pol 
end

outdir = joinpath(@__DIR__, "..", "data")
mkpath(outdir)   # ensure folder exists
outfile = joinpath(outdir, "Bfield_with_Traj_test.gif") 

t_start = input_parameters["min_time"] #0.0 * yr 
t_end = input_parameters["max_time"] # 1.0 * yr 
B_time_offset = input_parameters["B_field_time_offset"]

#times_number = Int( floor(27.0 / 5.0) * floor(365.0 / 27.0) * (t_end - t_start) / yr )  
times_number = 50 
times = range(t_start, t_end; length=times_number)  



gr()

n_dense = 500
vis_mode = "xz"
xspr = (0.1,105) 
yspr = xspr[2] - xspr[1]

anim = @animate for (i, t) in enumerate(times)

    @info "Rendering frame $i / $(length(times))  |  t = $(round(t/yr, digits=3)) yr"

    x, y, Bvals = B_data(
        input_parameters;
        radius = 50.0, #rad_val,
        timestamp   = t,
        n_x    = n_dense,
        n_t    = n_dense,
        n_z    = n_dense,
        visual_mode = vis_mode, 
        xspan_range = xspr, 
        yspan_range = yspr
    )

    if vis_mode == "xz"
        divisor = 1.0 

        title_string = "Time: $(round(t / yr, digits=2)) yr - " *
                        "$(input_parameters["particle_type"]) of $(input_parameters["particle_size"]) nm " * 
                        "with $(input_parameters["charging_type"]) charging"
        xlabel, ylabel = ("x [AU]", "z [AU]")
    else
        divisor = yr
        title_string = "Radius: $(rad_val) AU"
        xlabel, ylabel = ("t [yr]", "z [AU]")
    end

    plt = heatmap(
        x ./ divisor, y, Bvals', 
        margin = 8Plots.mm, size = (700, 800), title = title_string,
        xlabel = xlabel, ylabel = ylabel, colorbar = false, c=:balance )

    idx_time = argmin(abs.(sol_time .- t))

    plot!(plt, sol_x[1:idx_time] ./ AU, sol_z[1:idx_time] ./ AU,
        lw = 4, color = :red, label = "")

    scatter!(plt,
        [sol_x[idx_time]] ./ AU,
        [sol_z[idx_time]] ./ AU,
        color = :red,
        markersize = 3,
        label = "")
end

gif(anim, outfile; fps=2)