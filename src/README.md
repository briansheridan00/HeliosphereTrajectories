# Source Code Description

Calculating the trajectory of a charged interstellar dust particle requires the accurate modelling of the dust grain properties, the environment it passes through, and the forces acting on the particle. These forces are the gravitational force, solar radiation pressure force, and the electromagnetic Lorentz force. For particles far from the sun and interacting at the heliospheric boundary, the gravitational and solar radiation pressure forces are negligible, allowing us to focus on the Lorentz force. 

The magnetic field around the sun originates from the movement of conducting material within the sun. The creates a magnetic north and south pole with magnetic field lines stretching around the sun. The sun also emits a solar wind, which is a stream of charged particles. We assume that the magnetic field lines get frozen into this solar wind, therefore transporting the magnetic field outward into the solar system. There exists a boundary layer between the magnetic field lines directed toward or away from the sun since we have a south and north magnetic pole, defining positive and negative polarities and a resulting heliospheric current sheet separating these polarites. 

There are various ways of modelling these two phenomena for the solar wind and the magnetic field. Simplified approaches may be employed, such as a radial solar wind velocity and constant magnetic field values across the heliospheric current sheet boundary region. More nuanced models may incorporate a turning plasma velocity flow or a Parker spiral for the magnetic field, which are implemented within this code. The following table presents an overview and description of the relevant files to compute the trajectory in full. 

| File Name | Description |
| :--- | :--- | 
| `build_charging_data.jl` | A file to build a dictionary of charging times and equilibrium voltages for a range of different particle compositions, sizes, and regions of the heliosphere, as determined by previous work. **Note:** The code saves the dictionary as a JLD2 (HDF5) file at `/data/charging_dict.jld2`. |
| `charging.jl` | A collection of functions to interpolate between the saved values within the `charging_dict.jld2` file to determine particle properties, charging times, charge-to-mass ratios and equilibrium voltages. |
| `constants.jl` | A collection of relevant constants for use in the computations.  |
| `magnetic_field.jl` | Code to implement the magnetic field according to various prescribed models, including `simple`, `constant`, `solar_min`, `solar_max`, and `Parker`. Additional function to calculate the polarity for a grid of z-axis and time values for the Parker spiral mode, necessary for plotting. |
| `plasma_field.jl` | Implementation of various solar wind models, which differ in geometry. Available models are called `straight_straight`, `vertical_straight`, `vertical_turning`, and `turning_turning`. |
| `plotting.jl` | A collection of functions to plot the trajectory and varoius other physical information. Discussed in more detail below.  |
| `trajectory.jl` | Code to calculate the trjaectory of a particle given the input file parameters. There are three separate method used, namely `constant`, `instant`, and `continuous`. Discussed in more detail below. |
| `utils.jl` | A collection of useful computational functions and physics functions, which aid in the extraction of parameters for the trajectories as well as the calculation of physical quantities. Especially important is the `load_parameters()` funciton which should always be used to read out the input parameters file as a dictionary. |


## Plotting 

#### Main plotting function

The plotting capabilities of the code are twofold. The first is a overarching plotting function called `Plotter(sol, input; plot_B_fields=true, charges=false, n_density=300)`. This function produces a plot of the trajectory, and optional additional physical characteristics if required. Specifically, the main plot is a slice plot of the trajectory across the z-x plane (or optionally the y-x plane), with the particle path coloured by the speed at that point. The plasma field is denoted by quiver arrows showing the direction (magnitude ignored). A range of time steps for the trjaectory is included beside the particle path, allowing one to easily understand the behaviour over time. If the `plot_B_fields` argument is true, there is a triplet magnetic field heatmaps in the z-x plane at various timesteps. If the `charges` argument is true, then the charge to mass ratio is plotted for hte particle over time, which can change as it charges up or moves across boundaries. The `n_density` value determines the resolution of the magnetic field heatmaps. A wireframe diagram of the plot is available here. 

<pre> 
┌───────────────────────────────┬───────────┐
│                               │ Magnetic  │
│        Trajectory Plot        │ Field     │
│                               │ Triplet   │
│                               │           │
│                               │           │
├───────────────────────────────┬───────────┤
│   Charge-Mass Ratio vs Time   │ Magnetic  │
│                               │ Polarity  │
└───────────────────────────────┴───────────┘
</pre>

This function relies on a number of subfunctions to build up the constituent parts. These include `SlicePlot()` to downsample the trajectory and then plot the trjaecotry plot, where a number of other functions add the plasma field arrows, the termination shock and heliopause boundary lines, and the timesteps along the trajectory path. There are then additional functions to plot the heatmaps of the magnetic field at various timesteps, which draw on the code within `magnetic_field.jl`. Finally, the function `Bfield_over_time()` plots the magnetic field polarity as experienced by the particle as a function over time. 

#### Animated plots 

Animated plots are able to be created by creating a stack of images for a given particle trajectory. Once we compute the trajectory, we can compute the magnetic field for a range of timesteps and create heatmap of the magnetic field polarity and the trjaectory up to that point for each timestep. The stack of images is collected into an animated GIF. Care should be taken not to have a timestep between frames which is a significant fraction, or greater, than the solar rotation period, as this can create optical illusions within the GIF. The code for this is available within `main/animated_trajectories.jl`. 


## Trajectory  
The trajectory is computed using one of three different methods, either `constant`, `instant`, or `continuous`, where each method corresponds to the charging process of the particle. In the case of a constant charge to mass ratio, the $q/m$ value fo the particle never changes, and the trajectory is calculated using the standard Julia package `DifferentialEquations.jl`, as in every case. In contrast, instant charging refers to a sharp change in the $q/m$ value as the particle passes a boundary. During the particle trajectoy computation using `DifferentialEquations.jl`, we use callback functions to continuously monitor whether the pasticle crosses a boundar and in which direction. The charge to mass ratio is accordingly updated using the charging properties dictionary and utility functions from `utils.jl`. The third case refers to a continuous change in the particle $q/m$ along the trajectory. At each timestep the particle position is checked, the equilibrium $q/m$ for that region and particle characteristics found, and the an incremental change is added to the current $q/m$ in an exponentially changing manner. All these three charging methods, or lack thereof, allow flexibility of charging implementation and easy comparison between models. 

