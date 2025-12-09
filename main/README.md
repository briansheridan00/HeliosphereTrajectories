# Running Trajectories 

The input file contains the relevant simulation information to run the trajectory and plot the results. 

Note the following description of the input file parameters. 


| Parameter | Description |
| :--- | :--- | 
| `distance_HP` | Distance from the sun to the **Heliopause** (boundary of the heliosphere) in AU. |
| `distance_TS` | Distance from the sun to the **Termination Shock** in AU. |
| `B_mag_Heliopause` | Magnetic field strength ($B$) within the **Heliosheath** in Teslas (T). |
| `B_mag_ISM` | Magnetic field strength ($B$) in the **Interstellar Medium (ISM)** in Teslas (T). |
| `ISM_plasma_velocity` | Bulk plasma velocity of the **Interstellar Medium** flow in km/s. |
| `Heliosheath_plasma_velocity` | Bulk plasma velocity within the **Heliosheath** (region between TS and HP), in km/s. |
| `Termination_plasma_velocity` | Typical plasma velocity within the **Termination Shock** (Solar Wind), in km/s. |
| `distance_Turning` | Scale lenth of the turning of the plasma velocity field beyond the heliopause (in the ISM), in AU. |
| `distance_Approach` | Scale lenth of the turning of the plasma velocity field within the heliopause (in the heliosheath), in AU. |
| `mode` | Simulation mode; options are `"full"` or `"reduced"`, for all three forces (Gravity, Solar Radiation Pressure, Lorentz) or just Gravity, respectively. |
| `dist_measure` | Coordinate system geometry for measuring distance; options: `"spherical"` or `"flat"`. |
| `plasma_model` | The specific model used for the background plasma flow; options: `"straight_straight"`, `"vertical_straight"`, `"vertical_turning"`, `"turning_turning"`. | 
| `B_direction_Heliopause` | Direction of the magnetic field into the y-axis in the heliosheath; Options: `"+1"` or `"-1"`. |
| `B_model` | Magnetic field mode; Options: `"simple"`, `"constant"`, `"solar_min"`, `"solar_max"` or `"Parker"`. |
| `B_field_time_offset` | Time offset in years from the initial phase of the Parker spiral.|
| `x0_position` | **Initial position** $x_0$ of the particle in AU. |
| `y0_position` | **Initial position** $y_0$ of the particle in AU. |
| `z0_position` | **Initial position** $z_0$ of the particle in AU. |
| `alpha_angle` | **Initial angle** $\alpha$ defining the particle's initial velocity vector, in degrees, between x and z axes. |
| `beta_angle` | **Initial angle** $\beta$ defining the particle's initial velocity vector, in degrees, between x and y axes. |
| `min_time` | **Start time** for the trajectory integration in years. |
| `max_time` | **End time** or duration for the trajectory integration in years. |
| `dt` | The **time step** used for saving or outputting trajectory data in seconds. |
| `beta_value` | The ratio of solar radiation pressure force to gravitational force. |
| `q_over_m_value` | The **charge-to-mass ratio** ($q/m$) for the particle being tracked, in C/kg. |
| `qm_initial` | The **initial charge-to-mass ratio** ($q/m$) for the particle being tracked in the case of non-constant charge, in C/kg. |
| `particle_type` | The **particle material type**; Options: `"carbonaceous"` or `"silicate"` (with densitites 2.5 and 3.3 g/cm^3, respectively.) |
| `particle_size` | The **size of the particle** in nm, necessary to calculate equilibrium voltages, time-delays, and charge to mass ratios in various regions. |
| `charging_type` | The **method of particle charging**. Options; `"constant"`, `"instant"`, or `"continuous"`. Note for continuous charging, we currently implement exponential charging.  |
| `plane` | The plane which is plotted in two dimensions; Options: `"xz"`, `"xy"`, or `"both"` |
| `plot_sun` | The boolean value to choose whether to scatter plot the sun or not. |
| `plot_vectors` | The boolean value whether to plot the plasma velocity vectors or not. |
| `plot_magnetic_field` | The boolean value whether to plot the magnetic field direction or not. |
| `n_grid` | The number density of the grid of plasma velocity vector arrows. A suitable density is approximately 20. |
| `plot_color` | The value with which the trajectory is colored; Options: `"speed"` or `"time"`. |
| `save_fig` | A boolean option to save the figure in a data specific folder. |
| `annotate_time` | A boolean option to print time stamps in years nearby the trajectory. |