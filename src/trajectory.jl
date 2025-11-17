using DifferentialEquations 
using LinearAlgebra 
using TOML 
include(joinpath(@__DIR__, "constants.jl"))
include(joinpath(@__DIR__, "magnetic_field.jl"))
include(joinpath(@__DIR__, "plasma_field.jl")) 

#input_file_path = joinpath(@__DIR__, "..", "main", "input_values.toml") 
#input_dict = load_parameters(input_file_path) 


# Define the equation of motion 
# The ! indicates the function mutates its first argument (du) 
function EqMotion!(du, u, p, t, input) 

    """ 
    Define the equation of motion. 
    """
    mode = input["mode"]
    beta, q_over_m = input["p"] # p consists of the input parameters 

    # u = [position x, position y, position z, velocity x, velocity y, velocity z] 
    du[1:3] = u[4:6] 

    # Return updated acceleration equations 
    if mode == "reduced" # Just Gravity 
        rnorm = norm(u[1:3]) # Always have spherical rnorm for gravity. 
        du[4] = - (GM_Sun / rnorm^3) * u[1] 
        du[5] = - (GM_Sun / rnorm^3) * u[2] 
        du[6] = - (GM_Sun / rnorm^3) * u[3]  
        #du_vec = GravityForce(u) 
        #du[4], du[5], du[6] = du_vec 

    elseif mode == "full" # Gravity, Solar Radiation Pressure, and Lorentz Force 
        v_rel = u[4:6] - PlasmaVelocity(u; input=input) 
        B = B_field(u; input=input) 
        rnorm = norm(u[1:3]) # Always have spherical rnorm for gravity. 
        cross_v_B = cross(v_rel, B) 
        forces_prefactor = ((1 - beta) * GM_Sun / rnorm^3)
        du[4] = - forces_prefactor * u[1] + q_over_m * cross_v_B[1]
        du[5] = - forces_prefactor * u[2] + q_over_m * cross_v_B[2]
        du[6] = - forces_prefactor * u[3] + q_over_m * cross_v_B[3] 

    else # Error for invalid mode
        error("Invalid mode: '$mode'. Allowed modes are \"reduced\" or \"full\".")
    end 
end 


# Function to compute the trajectory solution - with keyword arguments. 
function ComputeTrajectory(input)
    # --- Extract parameters --- 
    params = input["p"] 

    # --- Compute angles in radians ---
    alpha_angle = deg2rad( input["alpha_angle"] )
    beta_angle  = deg2rad( input["beta_angle"] )

    # --- Initial velocity components ---
    v_mag = input["ISM_plasma_velocity"]
    vx0 = v_mag * cos(beta_angle) * cos(alpha_angle)
    vy0 = v_mag * cos(beta_angle) * sin(alpha_angle)
    vz0 = v_mag * sin(beta_angle)

    # --- Assemble initial condition ---
    r0 = input["r0"] 
    u0 = [r0[1], r0[2], r0[3], vx0, vy0, vz0]

    # --- Define time range ---
    min_time = input["min_time"]
    max_time = input["max_time"]
    tspan = (min_time, max_time)

    # --- Define ODE problem ---
    prob = ODEProblem(
        (du, u, p, t) -> EqMotion!(du, u, p, t, input), 
        u0, tspan, params 
    )

    # --- Solve the ODE ---
    sol = solve(prob, Vern9(), adaptive=false, dt = input["dt"])

    # --- Display summary ---
    println("--- Solution Characteristics ---")
    println("Integration time: $(min_time/yr) → $(max_time/yr) years")
    println("Solution size: $(size(sol))")

    return sol
end

