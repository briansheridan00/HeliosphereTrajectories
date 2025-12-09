using Pkg 
using TOML 
using DifferentialEquations 
using LinearAlgebra 
using JLD2
using FileIO 

include(joinpath(@__DIR__, "..", "src", "constants.jl"))
include(joinpath(@__DIR__, "..", "src", "magnetic_field.jl"))
include(joinpath(@__DIR__, "..", "src", "plasma_field.jl")) 
include(joinpath(@__DIR__, "..", "src", "charging.jl"))  

# --- Define struct for parameters: efficient for ODE solver --- 
#=
struct InputParams
    mode::String 
    beta::Float64
    qm_constant::Float64
    qm_initial::Float64 
    particle_type::String 
    particle_size::Float64 
    charging_dictionary::Dict 
end 
=#


# Function to compute the trajectory solution - with keyword arguments. 
function ComputeTrajectory(input)
    # --- Extract dictionary of dust particle charging data --- 
    charging_data = load(joinpath(@__DIR__, "..", "data", "charging_dict.jld2"))
    charging_dict = charging_data["charging_dict"] 

    # --- Extract parameters and conditions --- 
    mode = input["mode"] 
    charging_type = input["charging_type"] 
    beta_val = input["beta_value"] 
    qm_constant = input["q_over_m_value"] 
    qm_initial = input["qm_initial"] 
    particle_type = input["particle_type"]
    particle_size = input["particle_size"] 
    density = particle_type == "carbonaceous" ? 2.5 : 3.3

    # --- Calculate the initial q/m value for the ISM (case of instantaneous charging - no boundary crossed)  
    v_eq_ism = voltage_at_size(charging_dict, particle_type, "VLISM", particle_size)  
    qm_ism = calculate_qm(v_eq_ism, particle_size, density) 

    # --- Compute angles in radians ---
    alpha_angle = deg2rad( input["alpha_angle"] )
    beta_angle  = deg2rad( input["beta_angle"] )

    # --- Initial position and velocity components ---
    r0 = input["r0"] 
    v_mag = input["ISM_plasma_velocity"]
    vx0 = v_mag * cos(beta_angle) * cos(alpha_angle)
    vy0 = v_mag * cos(beta_angle) * sin(alpha_angle)
    vz0 = v_mag * sin(beta_angle) 

    # --- Define time range ---
    min_time = input["min_time"]
    max_time = input["max_time"]
    tspan = (min_time, max_time)

    # --- Define the affect and callback functions for instantaneous charging. 
    cond_HP(u,t,integrator) = norm(u[1:3]) - input["distance_HP"]
    cond_TS(u,t,integrator) = norm(u[1:3]) - input["distance_TS"]  

    function affectInstant!(integrator)
        rnorm_val = norm(integrator.u[1:3])  
        region = (rnorm_val ≤ input["distance_TS"]) ? "TerminationShock" : (rnorm_val ≤ input["distance_HP"]) ? "Heliosheath" : "VLISM" 
        v_eq = voltage_at_size(charging_dict, particle_type, region, particle_size)  
        qm_new = calculate_qm(v_eq, particle_size, density) 
        integrator.p[2] = qm_new # Q/m value  
    end  

    # Empty object to save the charge to mass ratio values. 
    saved_qm = SavedValues(Float64, Float64)

    # Function to push the q/m value to storage. 
    function save_qm(u, t, integrator)
        integrator.p[2]   # q/m in parameters
    end

    # Callback pushing to storage. 
    save_cb = SavingCallback(save_qm, saved_qm) 

    # Define callbak set with value saving to storage. 
    cb = CallbackSet(
        ContinuousCallback(cond_HP, affectInstant!),
        ContinuousCallback(cond_TS, affectInstant!), 
        save_cb 
    ) 

    # --- Define ODE problem and solve --- 
    if charging_type == "constant"
        params = (beta_val, qm_constant) 
        #params = InputParams(mode, beta_val, qm_constant, qm_initial, particle_type, particle_size, charging_dict)
        u0 = [r0[1], r0[2], r0[3], vx0, vy0, vz0] 
        prob = ODEProblem( (du, u, p, t) -> EqMotionConstant!(du, u, p, t, input), u0, tspan, params )
        sol = solve(prob, Vern9(), adaptive=false, dt = input["dt"])
    elseif charging_type == "instant"
        params = [beta_val, qm_ism] #qm_initial] 
        #params = InputParams(mode, beta_val, qm_constant, qm_initial, particle_type, particle_size, charging_dict)
        u0 = [r0[1], r0[2], r0[3], vx0, vy0, vz0]
        prob = ODEProblem( (du, u, p, t) -> EqMotionInstant!(du, u, p, t, input), u0, tspan, params )
        sol = solve(prob, Vern9(), adaptive=false, dt = input["dt"]; callback = cb)
    elseif charging_type == "continuous"
        params = (beta_val, particle_type, particle_size)
        #params = InputParams(mode, beta_val, qm_constant, qm_initial, particle_type, particle_size, charging_dict)
        u0 = [r0[1], r0[2], r0[3], vx0, vy0, vz0, qm_initial]
        prob = ODEProblem( (du, u, p, t) -> EqMotionContinuous!(du, u, p, t, input, charging_dict), u0, tspan, params )
        sol = solve(prob, Vern9(), adaptive=false, dt = input["dt"])
    else 
        error("Charging type not recognised")
    end 

    # --- Display summary ---
    println("--- Solution Characteristics ---")
    println("Integration time: $(min_time/yr) → $(max_time/yr) years")
    println("Solution size: $(size(sol))")


    # Return solution and charge values. 
    if charging_type == "constant"
        qm_array = fill(qm_constant, length(sol.t))

    elseif charging_type == "instant"
        # Extract callback-saved times + qm values
        saved_times = saved_qm.t
        saved_vals  = saved_qm.saveval
        # Prepare full Q/m array for every solver time
        qm_array = similar(sol.t, Float64)
        current_qm = qm_initial
        k = 1  # index into saved_vals
        for i in eachindex(sol.t)
            # Check if callback fired at this time (with tolerance)
            if k ≤ length(saved_times) && isapprox(sol.t[i], saved_times[k]; atol=1e-8)
                current_qm = saved_vals[k]
                k += 1
            end
            qm_array[i] = current_qm
        end

    elseif charging_type == "continuous"
        qm_array = sol[7, :]
    end

    return sol, qm_array

end



# --- Define the Equations of Motion - the ! indicates the function mutates its first argument (du) ---

# Constant charge to mass ratio. 
function EqMotionConstant!(du, u, p, t, input) 
    """ 
    Define the equation of motion. 
        Arguments: 
            u = [x, y, z, vx, vy, vz] 
            p = (beta, q_over_m) 
            t = (minimum_time, maximum_time) 
    """
    mode = input["mode"] 
    beta, q_over_m = p #[p.mode, p.beta, p.qm_constant] 

    # --- Positions and velocities contained within u. ---
    du[1:3] = u[4:6] 

    # --- Return updated acceleration equations ---
    if mode == "reduced" # Just Gravity 
        rnorm = norm(u[1:3]) # Always have spherical rnorm for gravity. 
        du[4:6] = - (GM_Sun / rnorm^3) .* u[1:3]  

    elseif mode == "full" # Gravity, Solar Radiation Pressure, and Lorentz Force 
        v_rel = u[4:6] - PlasmaVelocity(u, input) 
        B = B_field(u, input; t=t) 
        rnorm = norm(u[1:3]) # Always have spherical rnorm for gravity. 
        a_Gravity_SRP = - ((1 - beta) * GM_Sun / rnorm^3) .* u[1:3] 
        a_Lorentz = q_over_m .* cross(v_rel, B)  
        du[4:6] = a_Gravity_SRP + a_Lorentz 

    else # Error for invalid mode
        error("Invalid mode: '$mode'. Allowed modes are \"reduced\" or \"full\".")
    end 
end 


# Instantaneous change in charge to mass ratio upon crossing boundary. 
function EqMotionInstant!(du, u, p, t, input) 
    """ 
    Define the equation of motion. 
        Arguments: 
            u = [x, y, z, vx, vy, vz] 
            p = (beta, q_over_m) 
            t = (minimum_time, maximum_time) 
    """
    mode = input["mode"] 
    beta, q_over_m = p # Q/m is updated through callbacks. 

    du[1:3] = u[4:6] 

    if mode == "reduced"  
        rnorm = norm(u[1:3])  
        du[4:6] = - (GM_Sun / rnorm^3) .* u[1:3]  

    elseif mode == "full" 
        v_rel = u[4:6] - PlasmaVelocity(u, input) 
        B = B_field(u, input; t=t) 
        rnorm = norm(u[1:3]) 
        a_Gravity_SRP = - ((1 - beta) * GM_Sun / rnorm^3) .* u[1:3] 
        a_Lorentz = q_over_m .* cross(v_rel, B)  
        du[4:6] = a_Gravity_SRP + a_Lorentz 

    else 
        error("Invalid mode: '$mode'. Allowed modes are \"reduced\" or \"full\".")
    end 
end 


# Continuous change in charge to mass ratio. 
function EqMotionContinuous!(du, u, p, t, input, charge_d) 
    """ 
    Define the equation of motion. 
        Arguments: 
            u = [x, y, z, vx, vy, vz, qm_ratio] 
            p = (beta, type, size)
            t = (minimum_time, maximum_time) 
    """
    mode = input["mode"] 
    beta, type, size = p 
    density = type == "carbonaceous" ? 2.5 : 3.3 

    du[1:3] = u[4:6]  

    if mode == "reduced"  
        rnorm = norm(u[1:3])  
        du[4:6] = - (GM_Sun / rnorm^3) .* u[1:3] 
    elseif mode == "full" 
        v_rel = u[4:6] - PlasmaVelocity(u, input) 
        B = B_field(u, input; t=t) 
        rnorm = norm(u[1:3]) 
        a_Gravity_SRP = - ((1 - beta) * GM_Sun / rnorm^3) .* u[1:3] 
        a_Lorentz = u[7] .* cross(v_rel, B)  
        du[4:6] = a_Gravity_SRP + a_Lorentz  
        if rnorm <= input["distance_TS"] 
            region = "TerminationShock"
        elseif rnorm > input["distance_TS"] && rnorm <= input["distance_HP"]  
            region = "Heliosheath"
        else 
            region = "VLISM" 
        end  

        v_eq = voltage_at_size(charge_d, type, region, size)  
        tau = tau_value(charge_d, type, region, size) 
        q_eq = calculate_qm(v_eq, size, density) 
        du[7] = qm_dot_over_q(q_eq, u[7], tau) 
    else  
        error("Invalid mode: '$mode'. Allowed modes are \"reduced\" or \"full\".")
    end 
end 



