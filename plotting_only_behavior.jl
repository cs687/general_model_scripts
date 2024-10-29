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
include("../Codebase/behaviour_grid_loading.jl")
include("../Codebase/behaviour_grid_plotting.jl")

### Behaviour Map Parameters ###
S_grid = 10 .^(range(-1,stop=2,length=300))
D_grid = 10 .^(range(-1,stop=2,length=300))
τ_grid = [0.1,0.15,0.20,0.30,0.50,0.75,1.0,1.5,2.0,3.0,5.0,7.50,10.0,15.0,20.0,30.0,50.0,75.0,100.0]
v0_grid = [0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.20]
n_grid = [2.0,3.0,4.0]
η_grid = [0.001,0.002,0.005,0.01,0.02,0.05,0.1]

target_folder = "BehaviourMap"
dataset = DataSet(target_folder,S_grid,D_grid,τ_grid,v0_grid,n_grid,η_grid);


### Plot Settings ###

# Default values.
gr()
default(fmt = :png, framestyle=:box, grid=false, xguidefontsize=16, yguidefontsize=16, titlefontsize=18);

# Sets colour scales.
red_scale = make_red_scale(4);
bluegreen_scale = make_bluegreen_scale(4);
blue_scale= make_blue_scale(3);


####################
### Make Figure ####
####################

# Make behaviour grid plot.

#for tau in [0.1,0.2,], j in 1:8
params = [50,0.01,2,0.05]; #(τ, v0, n, η) might work: params = [0.1,0.02,3,0.1]
bg = BehaviourGrid(params,dataset); 
plot_behaviour_grid(bg,start_s_slice=100,idx_axis=false,xguide="D",yguide="S")
