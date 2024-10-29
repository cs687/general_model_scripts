####################
### Preparations ###
####################

### Generates Occurrence Counts ####

### Fetch Functions ###
include("../Codebase/basic_analysis.jl")
include("../Codebase/file_management.jl")
include("../Codebase/models.jl")
include("../Codebase/model_simulation.jl")
include("../Codebase/basic_plotting.jl")
include("../Codebase/utility.jl")
include("../Codebase/Functions_specific/gillespie_behaviour_grids.jl")

### Plot Settings ###

# Default values.
pythonplot();
default(fmt = :png, framestyle=:box, grid=false, xguidefontsize=16, yguidefontsize=16, titlefontsize=18);

# Sets colour scales.
red_scale = make_red_scale(4);
bluegreen_scale = make_bluegreen_scale(4);
blue_scale= make_blue_scale(3);


################
### Analysis ###
################

### Prepare Classification ###
begin
    v0,v,I,K,n,d,τ = [0.1 ,2.0, 1.0, 15.0, 3, 0.01, 100.0]
                                           #[v0, v, I, S,     D,     K, n, d, τ];
    parameters_no_activation =              [v0, v, I, 0.1 ,  150.0, K, n, d, τ];
    parameters_stochastic_pulsing =         [v0, v, I, 20.0,  150.0, K, n, d, τ];
    parameters_oscillation =                [v0, v, I, 100.0, 150.0, K, n, d, τ];
    parameters_stochastic_anti_pulsing =    [v0, v, I, 210.0, 150.0, K, n, d, τ];
    parameters_heterogeneous_activation =   [v0, v, I, 0.25,  0.01,  K, n, d, τ];
    parameters_homogeneous_activation =     [v0, v, I, 1.0,   0.01,  K, n, d, τ];
    parameters_stable_bistability =         [v0, v, I, 0.15,  0.01,  K, n, d, τ];

    #parameters_stochastic_switching =      [v0, 1.2, I, 0.02,  0.01,  1.0, n, d, τ];
end

S_grid = 10 .^(range(-2,stop=2,length=300))
D_grid = 10 .^(range(-2,stop=2,length=300))
param_grid = [[v0, v, I, S, D, K, n, d, τ] for S in S_grid, D in D_grid]
parameters = [v0, v, I, K, n, d, τ]

prod(size(param_grid))*5/3600


### Make and Plot Classification

@time grid_classification = classify_sd_space_fast(S_grid,D_grid,parameters)
plot_behaviour_grid(grid_classification, S_grid, D_grid)

#@time classification_data = ClassificationData.(param_grid)
#behaviour_grid = classify_classification_data.(classification_data)
#plot_behaviour_grid(behaviour_grid, S_grid, D_grid)

get_behaviour_grid(grid_classification)
using Serialization
serialize("../Data/Simulated_data/GillespieMap/gillespie_grid_classification.jls", grid_classification)
grid_classification = deserialize("../Data/Simulated_data/GillespieMap/gillespie_grid_classification.jls")
