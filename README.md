![HelioPaths Logo](HelioPaths-Full.png "Charting dust trajectories through the solar system")

# Heliosphere Trajectories 

Julia code to investigate the trajectories of interstellar dust particles through the solar system. Particular focus on the influence of the heliopause, heliosheath, and termination shock on the scattering and filtering of dust particles. 


Work carried out within the Astrodust group at ETH Zürich.  


### Requirements
- Julia 
- TOML.jl 
- DifferentialEquations.jl
- Interpolations.jl
- Plots.jl 
- JLD2.jl
- FileIO.jl  


### Installation

Clone the repository and instantiate the environment:

```bash
git clone https://github.com/briansheridan00/HeliosphereTrajectories.git
cd HeliosphereTrajectories
julia --project=.
```

Inside the Julia REPL, activate and instantiate the environment: 

```julia
julia> ] activate .
julia> ] instantiate 
```
