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
            B_vector = [0.0, Float64(input["B_direction_Heliopause"]) * input["B_mag_Heliopause"], 0.0]
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
        solar_period = 27.0 * 86400.0 # Seconds in a 27 day rotation period.  
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
        # Define the parameters  
        T_sun = 27.0 * 86400.0      # Rotation period of the sun (seconds). 
        Ω_sun = (2.0 * pi) / T_sun  # Angular velocity of solar rotation. 
        T_HCS = 11.0 * yr           # Period from one solar cycle to the next (seconds). 
        v_SW_TS = input["Termination_plasma_velocity"]
        v_SW_HP = input["Heliosheath_plasma_velocity"]
        v_SW_ISM = input["ISM_plasma_velocity"]
        r0 = 10.0 * R_sun           # Base for Parker spiral. 
        B0 = 2300.0e-9              # Base magnetic field strength. 
        dist_TS = input["distance_TS"]
        dist_HP = input["distance_HP"] 
        B_time_offset = input["B_field_time_offset"] 

        # Adjust time value 
        time = B_time_offset + t 

        # Trigonometric functions to define the longitude and latitude of the position vector. 
        r = sqrt(dot(u[1:3], u[1:3]))
        rρ = sqrt(dot(u[1:2], u[1:2])) 

        if r < r0
            B_vector = [0.0, 0.0, 0.0]
            return B_vector 
        end 

        sinϕ = u[2] / rρ
        cosϕ = u[1] / rρ 
        sinθ = u[3] / r
        cosθ = rρ / r  

        # Determine HCS and solar rotation axes wrt to xy plane and x-axis. 
        # ϵ_HCS: Angle of HCS with respect to xy-plane ±(0°,360°)
        # α_sun: Angle of solar rotation with respect to x-axis ±(0°,360°) 
        if r <= dist_TS 
            ϵ_HCS = (180.0 / T_HCS * (time - ((r - r0) / v_SW_TS))) % 360.0 
            α_sun = (360.0 / T_sun * (time - ((r - r0) / v_SW_TS))) % 360.0 
        elseif r > dist_TS && r <= dist_HP
            ϵ_HCS = (180.0 / T_HCS * (time - (((dist_TS - r0) / v_SW_TS) + ((r - dist_TS) / v_SW_HP)))) % 360.0 
            α_sun = (360.0 / T_sun * (time - (((dist_TS - r0) / v_SW_TS) + ((r - dist_TS) / v_SW_HP)))) % 360.0 
        else 
            B_vector = [0.0, input["B_mag_ISM"], 0.0] 
            return B_vector 
        end 

        vec_sun = angle2vector((α_sun + 90.0) % 360.0, 0.0)
        ϵ_HCS += ϵ_HCS < 0.0 ? 360.0 : 0.0 

        if 0.0 <= ϵ_HCS <= 90.0
            vec_HCS = angle2vector(α_sun, ϵ_HCS)
        elseif 90.0 < ϵ_HCS <= 180.0
            vec_HCS = angle2vector(180.0 + α_sun, 180.0 - ϵ_HCS)
        elseif 180.0 < ϵ_HCS <= 270.0
            vec_HCS = angle2vector(180.0 + α_sun, 180.0 - ϵ_HCS)
        elseif 270.0 < ϵ_HCS <= 360.0
            vec_HCS = angle2vector(α_sun, ϵ_HCS - 360.0)
        end

        vec_normal = cross(vec_HCS, vec_sun)
        angle2HCS = VectorAngle(vec_normal, u[1:3]) 

        if angle2HCS < 90.0       # We are above the HCS plane
            p = 1.0
        elseif angle2HCS > 90.0   # We are below the HCS plane
            p = - 1.0
        else                      # We are on the HCS plane 
            p = 0.0 
        end 

        Br = p * B0 * (r0 / r)^2 
        Bphi = p * B0 * Ω_sun / v_SW_HP * (r - r0) * (r0 / r)^2 * cosθ 

        B_vector = [Br * cosθ * cosϕ + Bphi * sinϕ,
                    Br * cosθ * sinϕ - Bphi * cosϕ, 
                    Br * sinθ] 
    end 
    
    return B_vector  
end 
