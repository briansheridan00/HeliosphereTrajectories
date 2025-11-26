using TOML 
using LinearAlgebra

include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "constants.jl"))


function PlasmaVelocity(u, input; t=0.0) 
    """
    Define the plasma velocty vector in each region. The magnetic 
        field is frozen into the plasma velocity field.  
    """  
    dist_measure                 = input["dist_measure"]
    plasma_model                 = input["plasma_model"]
    distance_HP                  = input["distance_HP"]
    distance_TS                  = input["distance_TS"]
    Heliosheath_plasma_velocity  = input["Heliosheath_plasma_velocity"]
    Termination_plasma_velocity  = input["Termination_plasma_velocity"] 
    ISM_plasma_velocity          = input["ISM_plasma_velocity"]
    turning_distance             = input["distance_Turning"]
    distance_Approach            = input["distance_Approach"] 

    rnorm = rnorm_func(u, input) 

    if plasma_model == "straight_straight" 
        if rnorm <= distance_HP && rnorm > distance_TS
            if dist_measure == "flat"
                Plasma_vec = [Heliosheath_plasma_velocity, 0.0, 0.0]  
            elseif dist_measure == "spherical"
                Plasma_vec = Heliosheath_plasma_velocity .*  u[1:3] / norm(u[1:3]) 
            else  
                Plasma_vec = [Heliosheath_plasma_velocity, 0.0, 0.0]
            end 
        elseif rnorm <= distance_TS  
            if dist_measure == "flat"
                Plasma_vec = [Termination_plasma_velocity, 0.0, 0.0]  
            elseif dist_measure == "spherical"
                Plasma_vec = Termination_plasma_velocity .*  u[1:3] / norm(u[1:3]) 
            else  
                Plasma_vec = [Termination_plasma_velocity, 0.0, 0.0]
            end 
        else
            Plasma_vec = [-ISM_plasma_velocity, 0.0, 0.0]            
        end
    

    elseif plasma_model == "vertical_straight"
        if rnorm <= distance_HP && rnorm > distance_TS 
            if dist_measure == "flat"
                Plasma_vec = Heliosheath_plasma_velocity .* [0.0, 0.0, u[3] >= 0 ? 1.0 : -1.0]
            elseif dist_measure == "spherical"
                r_vec = u[1:3]
                x, y, z = r_vec[1], r_vec[2], r_vec[3]  
                rho = sqrt(x^2 + y^2)
                if rho < 1e-12
                    theta_hat = [0.0, 0.0, 1.0]  # Particle is on z axis, fallback direction
                    theta_hat .= z > 0 ? -theta_hat : theta_hat
                else
                    theta_hat = [x*z/(rnorm*rho), y*z/(rnorm*rho), -rho/rnorm] # Polar unit vector in Cartesian coordinates 
                    theta_hat .= z > 0 ? -theta_hat : theta_hat                  
                end
                Plasma_vec = Heliosheath_plasma_velocity .* theta_hat
            else  
                Plasma_vec =  Heliosheath_plasma_velocity .* [0.0, 0.0, u[3] >= 0 ? 1.0 : -1.0]
            end 
        elseif rnorm <= distance_TS 
            if dist_measure == "flat"
                Plasma_vec = [Termination_plasma_velocity, 0.0, 0.0]
            elseif dist_measure == "spherical"
                Plasma_vec = Termination_plasma_velocity .*  u[1:3] / norm(u[1:3]) 
            else  
                Plasma_vec = [Termination_plasma_velocity, 0.0, 0.0] 
            end 
        else
            Plasma_vec = [-ISM_plasma_velocity, 0.0, 0.0]
        end    


    elseif plasma_model == "vertical_turning" 
        if rnorm <= distance_HP && rnorm > distance_TS 
            if dist_measure == "flat"
                Plasma_vec =  Heliosheath_plasma_velocity .* [0.0, 0.0, u[3] >= 0 ? 1.0 : -1.0]
            elseif dist_measure == "spherical"
                r_vec = u[1:3]
                x, y, z = r_vec[1], r_vec[2], r_vec[3] 
                rho = sqrt(x^2 + y^2)
                if rho < 1e-12
                    theta_hat = [0.0, 0.0, 1.0]  # Particle is on z axis, fallback direction
                    theta_hat .= z > 0 ? -theta_hat : theta_hat
                else
                    theta_hat = [x*z/(rnorm*rho), y*z/(rnorm*rho), -rho/rnorm] # Polar unit vector in Cartesian coordinates 
                    theta_hat .= z > 0 ? -theta_hat : theta_hat 
                end
                Plasma_vec = Heliosheath_plasma_velocity .* theta_hat
            else  
                Plasma_vec =  Heliosheath_plasma_velocity .* [0.0, 0.0, u[3] >= 0 ? 1.0 : -1.0]
            end 
        elseif rnorm <= distance_TS 
            if dist_measure == "flat"
                Plasma_vec = [Termination_plasma_velocity, 0.0, 0.0]
            elseif dist_measure == "spherical"
                Plasma_vec = Termination_plasma_velocity .*  u[1:3] / norm(u[1:3]) 
            else  
                Plasma_vec = [Termination_plasma_velocity, 0.0, 0.0] 
            end 
        elseif rnorm > distance_HP + turning_distance 
            Plasma_vec = [-ISM_plasma_velocity, 0.0, 0.0]
        else 
            x_temp = 1.0 - ( (rnorm - distance_HP) / turning_distance)
            turning_angle = x_temp * pi / 2.0 
            Plasma_vec_x = - ISM_plasma_velocity * cos(turning_angle)
            Plasma_vec_y = 0.0
            Plasma_vec_z = ISM_plasma_velocity * sin(turning_angle) * (u[3] < 0 ? -1.0 : 1.0)
            Plasma_vec = [Plasma_vec_x, Plasma_vec_y, Plasma_vec_z] 
        end             


    elseif plasma_model == "turning_turning"  
        if rnorm <= distance_HP && rnorm > distance_TS 
            approach_distance = distance_HP - rnorm

            if dist_measure == "flat"
                if approach_distance < distance_Approach 
                    x_temp = 1.0 - ( approach_distance / distance_Approach ) 
                    turning_angle = x_temp * pi / 2.0 
                    Plasma_vec_x = Heliosheath_plasma_velocity * cos(turning_angle)
                    Plasma_vec_y = 0.0 
                    Plasma_vec_z = Heliosheath_plasma_velocity * sin(turning_angle) * (u[3] < 0 ? -1.0 : 1.0)
                    Plasma_vec = [Plasma_vec_x, Plasma_vec_y, Plasma_vec_z] 
                else 
                    Plasma_vec = [Heliosheath_plasma_velocity, 0.0, 0.0] 
                end 

            elseif dist_measure == "spherical" 
                if approach_distance < distance_Approach  
                    r_vec = u[1:3]
                    x, y, z = r_vec[1], r_vec[2], r_vec[3] 
                    rho = sqrt(x^2 + y^2)
                    if rho < 1e-12 
                        theta_hat = [0.0, 0.0, 1.0] #.* sign(z + 1e-12) #polar dir 
                        theta_hat .= z < 0 ? -theta_hat : theta_hat 
                    else 
                        theta_hat = [x*z/(rnorm*rho), y*z/(rnorm*rho), -rho/rnorm] #polar in cartesian  
                        theta_hat .= z > 0 ? -theta_hat : theta_hat 
                    end 
                    x_temp = 1.0 - (approach_distance / distance_Approach)
                    turning_angle = x_temp * pi / 2.0 
                    r_hat = r_vec / rnorm 
                    Plasma_vec = Heliosheath_plasma_velocity * ( cos(turning_angle) .* r_hat + sin(turning_angle) .* theta_hat) 
                else 
                    Plasma_vec = Heliosheath_plasma_velocity .*  u[1:3] / norm(u[1:3]) 
                end 
            end 
        
        elseif rnorm <= distance_TS 
            if dist_measure == "flat"
                Plasma_vec = [Termination_plasma_velocity, 0.0, 0.0]
            elseif dist_measure == "spherical"
                Plasma_vec = Termination_plasma_velocity .*  u[1:3] / norm(u[1:3]) 
            else  
                Plasma_vec = [Termination_plasma_velocity, 0.0, 0.0] 
            end 

        elseif rnorm > distance_HP + turning_distance 
            Plasma_vec = [-ISM_plasma_velocity, 0.0, 0.0]
        
        else 
            x_temp = 1.0 - ( (rnorm - distance_HP) / turning_distance)
            turning_angle = x_temp * pi / 2.0 
            Plasma_vec_x = - ISM_plasma_velocity * cos(turning_angle)
            Plasma_vec_y = 0.0
            Plasma_vec_z = ISM_plasma_velocity * sin(turning_angle) * (u[3] < 0 ? -1.0 : 1.0)
            Plasma_vec = [Plasma_vec_x, Plasma_vec_y, Plasma_vec_z] 
        end 


    else 
        error("Plasma model not recognised")
    end 

    return Plasma_vec 
end 