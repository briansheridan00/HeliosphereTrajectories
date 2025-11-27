using Pkg 
using TOML  
using Statistics 
using Random  

include( joinpath(@__DIR__, "..", "src", "constants.jl") ) 
include( joinpath(@__DIR__, "..", "src", "utils.jl") ) 
include( joinpath(@__DIR__, "..", "src", "magnetic_field.jl") ) 
include( joinpath(@__DIR__, "..", "src", "plasma_field.jl") ) 
include( joinpath(@__DIR__, "..", "src", "plotting.jl") ) 
include( joinpath(@__DIR__, "..", "src", "trajectory.jl") ) 

input_file_path = joinpath(@__DIR__, "..", "main", "input_values.toml") 
input_dict = load_parameters(input_file_path) 

# Compute trajectory. 
sol = ComputeTrajectory( input_dict )

# Plot the trajectory by executing the function.  
Plotter(sol, input_dict; plot_B_fields=true)

