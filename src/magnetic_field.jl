using TOML 
using LinearAlgebra
include( joinpath(@__DIR__, "utils.jl") ) 
include( joinpath(@__DIR__, "constants.jl") ) 



function B_field(u, input; t=0.0) 
    """
    Define the magnetic field values in each region. 
    """ 
    rnorm = rnorm_func(u, input)
    dist_measure = input["dist_measure"] 
    magnetic_model = input["B_model"]

    if magnetic_model == "simple"  
        if rnorm <= input["distance_HP"]
            B_vector = [0.0,  Float64(input["B_direction_Heliopause"]) * input["B_mag_Heliopause"], 0.0]
        else 
            B_vector = [0.0, input["B_mag_ISM"], 0.0] 
        end 
    
    elseif magnetic_model == "constant"
        if dist_measure == "spherical" 
            if rnorm <= input["distance_HP"] 
                phi = atan( u[2], u[1] ) # This is arctan2(y,x). 
                #phi = ifelse(abs(phi) < 1e-16, 0.0, phi) # Snap to zero if very small. 
                B_mag = input["B_mag_Heliopause"]
                Bx = B_mag * ( -sin(phi) )
                By = B_mag * ( cos(phi) )
                Bz = 0.0 
                B_vector = [Bx, By, Bz] .* Float64(input["B_direction_Heliopause"]) 
            else
                B_vector = [0.0, input["B_mag_ISM"], 0.0] 
            end 
        else # flat or otherwise 
            if rnorm <= input["distance_HP"]  
                B_vector = [0.0,  Float64(input["B_direction_Heliopause"]) * input["B_mag_Heliopause"], 0.0]
            else
                B_vector = [0.0, input["B_mag_ISM"], 0.0] 
            end 
        end 


    elseif magnetic_model == "solar_min" 
        if dist_measure == "spherical" 
            if rnorm <= input["distance_HP"] 
                phase_val = Float64(input["B_direction_Heliopause"]) 
                phi = atan( u[2], u[1] ) #arctan2(y,x) 
                B_mag = input["B_mag_Heliopause"]
                Bx = B_mag * ( -sin(phi) )
                By = B_mag * ( cos(phi) ) 
                Bz = 0.0 
                B_vector = [Bx, By, Bz] .* ( u[3] >= 0 ? phase_val : -phase_val ) 
            else
                B_vector = [0.0, input["B_mag_ISM"], 0.0] 
            end 
        else # flat or otherwise 
            if rnorm <= input["distance_HP"]  
                phase_val = Float64(input["B_direction_Heliopause"])
                B_vector = [0.0, input["B_mag_Heliopause"] , 0.0] .* ( u[3] >= 0 ? phase_val : -phase_val )
            else
                B_vector = [0.0, input["B_mag_ISM"], 0.0] 
            end 
        end 


    elseif magnetic_model == "solar_max" #error("Solar maximum model not implemented yet") 
        solar_period = 27.0 * 86400 # Seconds in a 27 day rotation period.  
        Termination_plasma_velocity = input["Termination_plasma_velocity"]         
        Heliosheath_plasma_velocity = input["Heliosheath_plasma_velocity"]

        if dist_measure == "spherical"
            if rnorm <= input["distance_HP"] && rnorm > input["distance_TS"]  
                condition = (t - (input["distance_TS"] / Termination_plasma_velocity + (rnorm - input["distance_TS"]) / Heliosheath_plasma_velocity )) 
                phi = atan( u[2], u[1] ) #arctan2(y,x) 
                B_mag = input["B_mag_Heliopause"]
                Bx = B_mag * ( -sin(phi) )
                By = B_mag * ( cos(phi) ) 
                Bz = 0.0
                B_vector = [Bx, By, Bz] .* ( abs(condition) % solar_period <= solar_period/2.0 ? 1.0 : -1.0 )

            elseif rnorm <= input["distance_TS"] 
                condition = (t - (rnorm / Termination_plasma_velocity))
                phi = atan( u[2], u[1] ) #arctan2(y,x) 
                B_mag = input["B_mag_Heliopause"]
                Bx = B_mag * ( -sin(phi) )
                By = B_mag * ( cos(phi) ) 
                Bz = 0.0 
                B_vector = [Bx, By, Bz] .* ( abs(condition) % solar_period <= solar_period/2.0 ? 1.0 : -1.0 )

            elseif rnorm > input["distance_HP"]
                B_vector = [0.0, input["B_mag_ISM"], 0.0]  
            end  

        else  
            if rnorm <= input["distance_HP"] && rnorm > input["distance_TS"]   
                condition = (t - (input["distance_TS"] / Termination_plasma_velocity + (rnorm - input["distance_TS"]) / Heliosheath_plasma_velocity )) 
                B_vector = [0.0, input["B_mag_Heliopause"] * ( abs(condition) % solar_period <= solar_period/2.0 ? 1.0 : -1.0 ) , 0.0] 

            elseif rnorm <= input["distance_TS"]  
                condition = (t - (rnorm / Termination_plasma_velocity))
                B_vector = [0.0, input["B_mag_Heliopause"] * ( abs(condition) % solar_period <= solar_period/2.0 ? 1.0 : -1.0 ) , 0.0]

            elseif rnorm > input["distance_HP"]
                B_vector = [0.0, input["B_mag_ISM"], 0.0] 
            end 
        end 


    elseif magnetic_model == "Parker"
        error("Parker model not implemented yet")
    end 
    
    return B_vector  
end 
