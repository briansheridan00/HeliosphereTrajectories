# Running Trajectories 

...

Note the following description of the input file parameters. 


| Parameter | Description |
| :--- | :--- | 
| `distance_HP` | Distance from the sun to the **Heliopause** (boundary of the heliosphere) in AU. |
| `distance_TS` | Distance from the sun to the **Termination Shock** in AU. |
| `B_mag_Heliopause` | Magnetic field strength ($B$) at the Heliopause boundary in Teslas (T). |
| `B_mag_ISM` | Magnetic field strength ($B$) in the **Interstellar Medium (ISM)** in Teslas (T). |
| `ISM_plasma_velocity` | Bulk plasma velocity of the **Interstellar Medium** flow in km/s. |
| `Heliosheath_plasma_velocity` | Bulk plasma velocity within the **Heliosheath** (region between TS and HP), in km/s. |
| `Inner_Solar_system_velocity` | Typical plasma velocity in the **Inner Solar System** (Solar Wind), in km/s. |
| `distance_Turning` | Scale lenth of the turning of the plasma velocity field, in AU. |
| `mode` | Simulation mode; options are `"full"` or `"reduced"`, for all three forces (Gravity, SOlar Radiation Pressure, Lorentz) or just Gravity. |
| `dist_measure` | Coordinate system for measuring distance; options: `"spherical"` or `"flat"`. |
| `plasma_model` | The specific model used for the background plasma flow; options: `"straight_straight"`, `"vertical_straight"`, `"vertical_turning"`, `"turning_turning"`. | 
| `x0_position` | **Initial position** $x_0$ of the particle in AU. |
| `y0_position` | **Initial position** $y_0$ of the particle in AU. |
| `z0_position` | **Initial position** $z_0$ of the particle in AU. |
| `alpha_angle` | **Initial angle** $\alpha$ defining the particle's initial velocity vector, in degrees, between x and z planes. |
| `beta_angle` | **Initial angle** $\beta$ defining the particle's initial velocity vector, in degrees, between x and y planes. |
| `min_time` | **Start time** for the trajectory integration in years. |
| `max_time` | **End time** or duration for the trajectory integration in years. |
| `dt` | The **time step** used for saving or outputting trajectory data in seconds. |
| `beta_value` | The ratio of solar radiation pressure force to gravitational force. |
| `q_over_m_value` | The **charge-to-mass ratio** ($q/m$) for the particle being tracked, in C/kg. |
| `plane` | The plane which is plotted din two dimensions; Options: `"xz"` or `"xy"` |
| `plot_sun` | The boolean value to choose whether to plot the sun or not. |
| `plot_vectors` | The boolean value whether to plot the plasma velocity vectors or not. |
| `n_grid` | The number density of the grid of plasma velocity vector arrows. Often choose 20. |
| `plot_color` | The value with which the trajectory is colored; Options: "speed" or "time". |
| `save_fig` | A boolean option to save the figure in a data specific folder. |