using Pkg
using TOML 
using LinearAlgebra
using Interpolations

include(joinpath(@__DIR__, "..", "src", "utils.jl"))
include(joinpath(@__DIR__, "..", "src", "constants.jl"))


# Calculate the Q/m value from particle and environment parameters. 
function calculate_qm(voltage, size, density) # Input units: [Volts, nm, g/cm^3] 
    e0 = 8.854188e-12        # Permittivity of free space. 
    U = voltage              # Voltage: in V. 
    rho = density * 1000     # Density: convert from g/cm^3 to kg/m^3. 
    rad = size * 1e-9        # Particle Radius: convert from nm to m. 
    qm = (3 * e0 * U) / (rho * rad^2) 
    return qm 
end 

# Function to get interpolated voltage for arbitrary sizes (single value or array)
function voltage_at_size(dict, material::String, region::String, size)
    # Extract the dust grain sizes. 
    charging_dict_sizes = dict["metadata"]["sizes"] 
    x = Float64.(charging_dict_sizes)

    # Extract discrete voltages
    volts = [dict[material][s][region]["Voltage"] for s in x ] 
    
    # Create linear interpolation function
    itp = LinearInterpolation(x, volts, extrapolation_bc=Line()) 
    
    # Broadcast interpolation over input (works for scalar or array)
    return itp.(size)
end

# Function to get interpolated charging time to 90 percent of equilibrium for arbitrary sizes (single value or array)
function charging_time_at_size(dict, material::String, region::String, size)
    # Extract the dust grain sizes. 
    charging_dict_sizes = dict["metadata"]["sizes"] 
    x = Float64.(charging_dict_sizes) 

    # extract the charging times for each size
    t90 = [dict[material][s][region]["ChargingTime"] for s in x]

    # linear interpolation
    itp = LinearInterpolation(x, t90, extrapolation_bc=Line())

    return itp.(size)
end

# Charge to mass ratio as a function of time. 
function qm_over_time(t, q_eq, q_0, tau)
    exp_term = 1 - exp(-t / tau)
    return q_0 + (q_eq - q_0) * exp_term 
end  

# Charge to mass ratio derivate as a function of charge. 
function qm_dot_over_q(q_eq, q_0, tau)
    return (q_eq - q_0) / tau 
end 

# Function to extract tau from the charging time to 90 percent. 
tau_from_t90(t90) = t90 / log(10) # Mathematically exact. 

# Fit for tau 
function tau_value(dict, material::String, region::String, size)
    time90 = charging_time_at_size(dict, material, region, size)
    return tau_from_t90(time90)
end