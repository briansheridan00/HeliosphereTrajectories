using TOML 
using LinearAlgebra
include(joinpath(@__DIR__, "constants.jl"))


function load_parameters(path::String)
    """ 
    Load a TOML parameter file and return a nested Dict.
    Distance-like variables are multiplied by `AU`.
    Time-like variables are multiplied by `yr`.
    """    
    # Load the raw TOML data (nested Dict)
    config = TOML.parsefile(path)

    config["distance_HP"] *= AU 
    config["distance_TS"] *= AU 
    config["distance_Turning"] *= AU  
    
    config["r0"] = [ config["x0_position"], config["y0_position"], config["z0_position"] ] .* AU  
    config["x0_position"] *= AU
    config["y0_position"] *= AU   
    config["z0_position"] *= AU 

    config["min_time"] *= yr 
    config["max_time"] *= yr 

    config["p"] = [config["beta_value"], config["q_over_m_value"]] 

    return config
end



### Function to calculate the norm of the position vector. 
function rnorm_func(u, input) 
    """ 
    Calculate the norm of the position vector from the sun if the geometry is flat or spherical. 
    """
    dist_measure = input["dist_measure"]

    if dist_measure == "flat" 
        return u[1] 
    elseif dist_measure == "spherical"
        return norm(u[1:3]) 
    else 
        @warn "Distance measure not recognised - defaulting to flat"
        return u[1] 
    end 
end 