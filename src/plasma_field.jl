using TOML 
using LinearAlgebra
include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "constants.jl"))

#input_file_path = joinpath(@__DIR__, "..", "main", "input_values.toml") 
#input_dict = load_parameters(input_file_path) 


function PlasmaVelocity(u, input) 
    """
    Define the plasma velocty vector in each region. The magnetic 
        field is frozen into the plasma velocity field.  
    """  
    dist_measure = input["dist_measure"]
    plasma_model = input["plasma_model"]
    distance_HP = input["distance_HP"]
    Heliosheath_plasma_velocity = input["Heliosheath_plasma_velocity"]
    ISM_plasma_velocity = input["ISM_plasma_velocity"]
    turning_distance = input["distance_Turning"]

    rnorm = rnorm_func(u, input) 

    if plasma_model == "straight_straight" 
        if rnorm <= distance_HP
            if dist_measure == "flat"
                Plasma_vec = [Heliosheath_plasma_velocity, 0.0, 0.0]  
            elseif dist_measure == "spherical"
                Plasma_vec = Heliosheath_plasma_velocity .*  u[1:3] / norm(u[1:3]) 
            else 
                @warn "Distance measure not recognised - plasma velocity falling back to flat" 
                Plasma_vec = [Heliosheath_plasma_velocity, 0.0, 0.0]
            end 
        #=
        elseif rnorm > distance_HP && rnorm <= Float64(distance_HP + turning_distance)
            x_temp = 1.0 - ( (rnorm - distance_HP) / turning_distance)
            turning_angle = x_temp * pi / 2.0 
            Plasma_vec = [- ISM_plasma_velocity * cos(turning_angle), 0, ISM_plasma_velocity * sin(turning_angle)]
        =# 
        else  #elseif rnorm > Float64( distance_HP + turning_distance ) 
            Plasma_vec = [-ISM_plasma_velocity, 0.0, 0.0]
        end
    

    elseif plasma_model == "vertical_straight"
        if rnorm <= distance_HP
            if dist_measure == "flat"
                Plasma_vec = Heliosheath_plasma_velocity .* [0.0, 0.0, u[3] >= 0 ? 1.0 : -1.0]
            elseif dist_measure == "spherical"
                r_vec = u[1:3]
                x, y, z = r_vec[1], r_vec[2], r_vec[3] #r_vec...
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
                @warn "Distance measure not recognised - plasma velocity falling back to flat" 
                Plasma_vec =  Heliosheath_plasma_velocity .* [0.0, 0.0, u[3] >= 0 ? 1.0 : -1.0]
            end                 
        else
            Plasma_vec = [-ISM_plasma_velocity, 0.0, 0.0]
        end    


    elseif plasma_model == "vertical_turning" 
        if rnorm <= distance_HP
            if dist_measure == "flat"
                Plasma_vec =  Heliosheath_plasma_velocity .* [0.0, 0.0, u[3] >= 0 ? 1.0 : -1.0]
            elseif dist_measure == "spherical"
                r_vec = u[1:3]
                x, y, z = r_vec[1], r_vec[2], r_vec[3] #r_vec...
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
                @warn "Distance measure not recognised - plasma velocity falling back to flat" 
                Plasma_vec =  Heliosheath_plasma_velocity .* [0.0, 0.0, u[3] >= 0 ? 1.0 : -1.0]
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
        error("Plasma model 'turning_turning' not implemented yet.") 
        
        
    else 
        error("Plasma model not recognised")
    end 


    return Plasma_vec 
end 