using TOML 
using LinearAlgebra
include(joinpath(@__DIR__, "constants.jl"))


# --- Trajectory functions --- # 
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



# --- Julia convenience functions --- #
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



function pretty_print_dict_old(d::Dict; indent::Int=0)
    prefix = " " ^ indent
    for (k, v) in d
        if isa(v, Dict)
            println("$prefix$k = {")
            pretty_print_dict(v; indent=indent+4)
            println("$prefix}")
        else
            println("$prefix$k = $v")
        end
    end
end
