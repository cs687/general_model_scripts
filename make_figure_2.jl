####################
### Preparations ###
####################

### Generates Occurrence Counts ####

### Fetch Functions ###
include("../Codebase/basic_analysis.jl")
include("../Codebase/file_management.jl")
include("../Codebase/utility.jl")
include("../Codebase/basic_plotting.jl")
include("../Codebase/models.jl")
include("../Codebase/model_simulation.jl")

### Plot Settings ###

# Default values.
# pythonplot();
default(fmt = :png, framestyle=:box, grid=false, xguidefontsize=16, yguidefontsize=16, titlefontsize=18);

# Sets colour scales.
red_scale = make_red_scale(4);
bluegreen_scale = make_bluegreen_scale(4);
blue_scale= make_blue_scale(3);


####################
### Make Figure ####
####################

lw = 3 # line width of plot.
t_pts = (-100.0, 0.0, 500.0) # (simulation starting point, simulation stress addition point, simulation end point)
ymax = 1.25; # maximum y value in simulation plot.

parameters_no_activation =              [1.0,  10.0, 10.0, 0.05,  2, 0.025]
parameters_stochastic_pulsing =         [2.8,  10.0, 10.0, 0.05,  2, 0.025]
parameters_oscillation =                [10.0, 10.0, 10.0, 0.05,  2, 0.025]
parameters_stochastic_anti_pulsing =    [15.5, 10.0, 10.0, 0.05,  3, 0.025]
parameters_heterogeneous_activation =   [2.25, 0.1,  10.0, 0.05,  2, 0.025]
parameters_homogeneous_activation =     [5.0,  0.1,  10.0, 0.05,  2, 0.025]
parameters_stochastic_switching =       [2.25, 1.2,  10.0, 0.05,  2, 0.025]
parameters_stable_bistability =         [2.1,  0.1,  10.0, 0.05,  2, 0.025]
parameters_single_pulse_response =      [2.75, 10.0, 10.0, 0.05,  2, 0.01]
parameters_intermediate_activation =    [8.0,  10.0, 0.1,  0.05,  2, 0.025]

# No activation.
# idxs is the idxs we plot (in [σ, A1, A2, A3], so `idx = [1]` is plotting sigma only).
plot_no_activation = plot_activation(parameters_no_activation, t_pts; idxs = [1], ymax=ymax, lw=lw)

# Stochastic pulsing.
plot_stochastic_pulsing = plot_activation(parameters_stochastic_pulsing, t_pts; idxs = [1], ymax=ymax, lw=lw)

# Oscillation.
plot_oscillation = plot_activation(parameters_oscillation, t_pts; idxs = [1], ymax=ymax, lw=lw)

# Stochastic anti-pulsing.
plot_stochastic_anti_pulsing = plot_activation(parameters_stochastic_anti_pulsing, t_pts; idxs = [1], ymax=ymax, lw=lw)

# Heterogeneous activation.
plot_heterogeneous_activation = plot_activations(parameters_heterogeneous_activation, t_pts, 3; ymax=ymax, lw=lw)

# Homogeneous activation.
plot_homogeneous_activation = plot_activation(parameters_homogeneous_activation, t_pts; idxs = [1], ymax=ymax, lw=lw)

# Stochastic switching.
plot_stochastic_switching = plot_activation(parameters_stochastic_switching, t_pts; idxs = [1], ymax=ymax, lw=lw)

# Stable bistability.
# Plot simulation simulation just a normal simulation plot.
# Here `(0.0, 500.0)` is simulation start and end point. `u0=[0.85, 0.85, 0.85, 0.85]` is which values we start simulation from.
plot_simulation(parameters_stable_bistability, (0.0, 500.0); u0=[0.85, 0.85, 0.85, 0.85], ymax=ymax, color=blue_scale[1], la=0.8, lw=lw)
plot_stable_bistability = plot_activation!(parameters_stable_bistability, t_pts; idxs = [1], ymax=ymax, color=blue_scale[2], la=0.8, lw=lw)

# Single pulse response.
plot_single_pulse_response = plot_activation(parameters_single_pulse_response, t_pts; idxs = [1], ymax=ymax, lw=lw)

# Intermediate activation.
plot_intermediate_activation = plot_activation(parameters_intermediate_activation, t_pts; idxs = [1], ymax=ymax, lw=lw)

# Make plot.
behaviours_figure = plot(plot_no_activation,plot_stochastic_pulsing,plot_oscillation,plot_stochastic_anti_pulsing,plot_heterogeneous_activation,plot_homogeneous_activation,plot_stochastic_switching,plot_stable_bistability,plot_single_pulse_response,plot_intermediate_activation,layout=(2,5),size=(3500,1000),left_margin=20mm,bottom_margin=25mm)    
make_figure(behaviours_figure, "Figure_2")
