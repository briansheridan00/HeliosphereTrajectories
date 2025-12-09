using Pkg
using TOML 
using JLD2 
using FileIO 

include( joinpath(@__DIR__, "..", "src", "constants.jl"))

# Toggle the following to save the dictionary. Note: Dictionary is only saved if file does not exist. 
save_dict = false  

# Look-up table of voltages and charging times. 
charging_dict = Dict() 

# Define the different data. 
charging_dict_materials = ["silicate", "carbonaceous"]
charging_dict_sizes = [20, 50, 100, 250, 500, 1000] # In nanometers 
charging_dict_regions = ["VLISM", "Heliosheath", "TerminationShock"]
charging_dict_properties = ["Voltage", "ChargingTime"] 

# Function to build the recursive dictionary. 
function build_recursive_dict(dimensions, default=nothing)
    if isempty(dimensions)
        return default
    end
    result = Dict()
    for key in first(dimensions)
        result[key] = build_recursive_dict(dimensions[2:end], default)
    end
    return result
end

# Define the structure of the dictionary. Outlines the nesting pattern. 
charging_dict_structure = [ 
    charging_dict_materials, 
    charging_dict_sizes, 
    charging_dict_regions, 
    charging_dict_properties 
    ]

# Build the data dictionary. 
charging_dict = build_recursive_dict(charging_dict_structure)  


"""
Data Information 
-> Structure: (x6 sizes, x3 regions) 
-> Sizes order: 20, 50, 100, 250, 500, 1000
-> Regions order: VLISM, Heliosheath, TerminationShock 
"""
silicate_voltages = [[1.37, 13.98, 21.12], [1.05, 12.37, 18.26], [0.78, 11.95, 16.16], 
                     [0.51, 11.74, 14.49], [0.38, 11.68, 13.89], [0.29, 11.66, 13.58] ] 

carbonaceous_voltages = [[2.35, 7.51, 18.8], [1.72, 6.8, 16.5], [1.03, 6.47, 14.59], 
                         [0.34, 6.22, 12.91], [0.05, 6.1, 12.27], [-0.09, 6.04, 11.95] ] 

silicate_times = [[4.3, 177.0, 1289.9], [2.0, 65.4, 591.0], [1.2, 32.0, 314.2], 
                  [0.7, 12.7, 128.2], [0.4, 6.3, 64.2], [0.2, 3.2, 32.1] ] .* day_value # convert days to seconds 

carbonaceous_times = [[8.5, 35.5, 974.7], [4.0, 13.5, 413.6], [2.6, 7.4, 200.8], 
                      [1.3, 3.3, 75.2], [0.7, 1.7, 36.3], [0.3, 0.9, 17.8] ] .* day_value # convert days to seconds 


# Function to add the data into the dictionary. 
function populate_data!(
    dict,
    material::String,
    voltages,
    times
)
    for (i, size) in enumerate(charging_dict_sizes)
        for (j, region) in enumerate(charging_dict_regions)
            dict[material][size][region]["Voltage"]      = voltages[i][j]
            dict[material][size][region]["ChargingTime"] = times[i][j]
        end
    end
end

populate_data!(charging_dict, "silicate", silicate_voltages, silicate_times)
populate_data!(charging_dict, "carbonaceous", carbonaceous_voltages, carbonaceous_times)

# Save the parameter data as metadata. 
charging_dict["metadata"] = Dict(
    "materials" => charging_dict_materials,
    "sizes"     => charging_dict_sizes,
    "regions"   => charging_dict_regions,
    "properties"=> charging_dict_properties
)

# If true and file does not exist, save as a JLD2 file - subset of HDF5. 
if save_dict == true
    save_path = joinpath(@__DIR__, "..", "data", "charging_dict.jld2")
    if !isfile(save_path)
        @save save_path charging_dict
        println("Saved charging_dict to: $save_path")
    else
        println("File already exists, skipping save: $save_path")
    end
end


"""
To load the dictionary, we use: 
> using JLD2, FileIO
> data = load(joinpath(@__DIR__, "..", "data", "charging_dict.jld2"))
> charging_dict = data["charging_dict"]
"""
