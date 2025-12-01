using TOML 
using LinearAlgebra
include(joinpath(@__DIR__, "constants.jl"))



### --- Physics functions --- ###

# Function to calculate the norm of the position vector. 
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

# Calculate the Q/m value from particle and environment parameters. 
function calculate_qm(voltage, radius, density) # Input units: [Volts, nm, g/cm^3] 
    e0 = 8.854188e-12        # Permittivity of free space. 
    U = voltage              # Voltage: in V. 
    rho = density * 1000     # Density: convert from g/cm^3 to kg/m^3. 
    rad = radius * 1e-9      # Radius: convert from nm to m. 
    qm = (3 * e0 * U) / (rho * rad^2) 
    return qm 
end 

# Angle between two vectors. 
function VectorAngle(u::AbstractVector, v::AbstractVector)
    cosinusθ = dot(u, v) / (norm(u) * norm(v))
    cosθ = clamp(cosinusθ, -1.0, 1.0)
    return acos(cosθ) * 180.0 / pi 
end

# Convert a pair of angles to a 3D unit vector. 
function angle2vector(phi_raw, theta_raw) 
    phi = phi_raw * pi / 180.0 
    theta = (90.0 - theta_raw) * pi / 180.0 
    return [ sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) ] 
end 

# Convert a vector direction to a pair of angles. 
function vector2angle(v)  
    phi = (sign(v[2]) * acos(v[1] / sqrt(v[1]^2 + v[2]^2))) / pi * 180.0
    theta = 90.0 - (acos(v[3] / sqrt(v[1]^2 + v[2]^2 + v[3]^2)) / pi * 180.0)
    return [ phi, theta ] 
end 


### --- Julia convenience functions --- ###

# Read a toml file and load the parameters as a Dict. Adjust parameter values if necessary. 
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
    config["distance_Approach"] *= AU 
    
    config["r0"] = [ config["x0_position"], config["y0_position"], config["z0_position"] ] .* AU  
    config["x0_position"] *= AU
    config["y0_position"] *= AU   
    config["z0_position"] *= AU 

    config["min_time"] *= yr 
    config["max_time"] *= yr 
    config["B_field_time_offset"] *= yr 

    config["p"] = [config["beta_value"], config["q_over_m_value"]] 

    return config
end

# Function to flatten nested dicts into key-value pairs with hierarchical keys
function flatten_dict(d::Dict, prefix::String="")
    rows = []
    for (k, v) in d
        full_key = isempty(prefix) ? string(k) : "$(prefix).$k"
        if isa(v, Dict)
            append!(rows, flatten_dict(v, full_key))
        else
            push!(rows, (full_key, v))
        end
    end
    return rows
end

# Function to print in a multi-column table
function pretty_print_table(d::Dict; ncols::Int=2)
    rows = flatten_dict(d)
    n = length(rows)
    nrows = ceil(Int, n / ncols)
    # Arrange items into table
    table_rows = ["" for _ in 1:nrows, _ in 1:ncols]
    for (i, (k, v)) in enumerate(rows)
        r = mod1(i, nrows)
        c = ceil(Int, i / nrows)
        table_rows[r, c] = rpad("$k = $v", 60)
    end
    # Print table
    for r in 1:nrows
        println(join(table_rows[r, :], "  "))
    end
end 

