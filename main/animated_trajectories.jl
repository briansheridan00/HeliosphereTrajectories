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

# Load the charging dictionary which interpolates between values for voltage and charging time. 
using JLD2, FileIO
data = load(joinpath(@__DIR__, "..", "data", "charging_dict.jld2"))
charging_dict = data["charging_dict"]

# Load the input file parameters using the utility function. 
input_file = joinpath(@__DIR__, "..", "main", "input_values.toml")
input_parameters = load_parameters(input_file) 

# Run the trajectory. 
res_traj, saved_charges = ComputeTrajectory(input_parameters); 

# Downsample the solution if there are too many points. 
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

# Define the path for the output file. 
outdir = joinpath(@__DIR__, "..", "data")
mkpath(outdir)   # ensure folder exists
outfile = joinpath(outdir, "Bfield_with_Traj_test.gif") 

# Define the starting time, ending time, and offset time for the magnetic field. 
t_start = input_parameters["min_time"] #0.0 * yr 
t_end = input_parameters["max_time"] # 1.0 * yr 
B_time_offset = input_parameters["B_field_time_offset"]

# Define the number and range of time steps. 
#times_number = Int( floor(27.0 / 5.0) * floor(365.0 / 27.0) * (t_end - t_start) / yr )  # Tuned to be every 5 days. 
times_number = 50 # Manually inputted amount. 
times = range(t_start, t_end; length=times_number)  



### --- Create the animation --- #
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

    # Plot the magnetic field solution 
    plt = heatmap(
        x ./ divisor, y, Bvals', 
        margin = 8Plots.mm, size = (700, 800), title = title_string,
        xlabel = xlabel, ylabel = ylabel, colorbar = false, c=:balance )

    # Extract the time of this frame and plot the trajectory up until this point. 
    idx_time = argmin(abs.(sol_time .- t))

    plot!(plt, sol_x[1:idx_time] ./ AU, sol_z[1:idx_time] ./ AU,
        lw = 4, color = :red, label = "")

    # Scatter plot the particle position at this time. 
    scatter!(plt,
        [sol_x[idx_time]] ./ AU,
        [sol_z[idx_time]] ./ AU,
        color = :red,
        markersize = 3,
        label = "")
end

# Save the gif in in the predefined destination. 
gif(anim, outfile; fps=2)
