using DifferentialEquations 
using LinearAlgebra 
using TOML 

include(joinpath(@__DIR__, "constants.jl"))
include(joinpath(@__DIR__, "magnetic_field.jl"))
include(joinpath(@__DIR__, "plasma_field.jl")) 


# Define the equation of motion - the ! indicates the function mutates its first argument (du). 
function EqMotion!(du, u, p, t, input) 
    """ 
    Define the equation of motion. 
        Arguments: 
            u = [x, y, z, vx, vy, vz] 
            p = [beta, q_over_m] 
            t = (minimum_time, maximum_time) 
    """
    mode = input["mode"]
    beta, q_over_m = input["p"] # p consists of the input parameters 

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

    # --- Assemble initial conditions ---
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

