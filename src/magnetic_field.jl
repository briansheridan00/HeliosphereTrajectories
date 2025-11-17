using TOML 
using LinearAlgebra
include( joinpath(@__DIR__, "utils.jl") ) 
include( joinpath(@__DIR__, "constants.jl") )

#input_file_path = joinpath(@__DIR__, "..", "main", "input_values.toml") 
#input_dict = load_parameters(input_file_path) 


function B_field(u, input) 
    """
    Define the magnetic field values in each region. 
    """ 

    rnorm = rnorm_func(u, input)

    if rnorm <= input["distance_HP"]  
       B_vector = [0, Float64(input["B_direction_Heliopause"]) * input["B_mag_Heliopause"], 0.0]
    else
        B_vector = [0, input["B_mag_ISM"], 0.0] 
    end 

    return B_vector  
end 
