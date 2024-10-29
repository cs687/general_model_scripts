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
params = [2.0,0.05,3,0.1];
bg = BehaviourGrid(params,dataset);
plot_behaviour_grid(bg,start_s_slice=100,idx_axis=false,xguide="D",yguide="S")

samples = [
    (0.2, 1.5),(0.2, 2.3),(0.2, 3.0),(0.2, 6.5),
    (1.5, 3.0),
    (3.5, 3.0),(3.5, 6.5),
    (10.0, 1.5),(10.0, 3.0),(10.0, 6.5),(10.0, 19.0),(10.0, 85.0),
    (60.0, 3.0),(60.0, 6.5),(60.0, 19.0),(60.0, 85.0)
]
sample_idxs = [(findfirst(D_grid .> pos[1]), findfirst(S_grid .> pos[2]) - 100) for pos in samples]
behaviour_plot = scatter!(sample_idxs,label="",color=:white,xlimit=(-Inf,Inf),ylimit=(-Inf,Inf))

# Make simulation plots.
make_p(ps, sample) = [sample[2]; sample[1]; ps]
sample_ps = [make_p(params, sample) for sample in samples]

lw = 3
t_pts = (-50.0, 0.0, 250.0)
ymax = 1.35;
plots = [plot_activation(p, t_pts; idxs = [1], ymax=ymax, lw=lw, xguide="", yguide="", saveat=0.2) for p in sample_ps]

# Handles special case: Heterogeneous activation.
i = 3
plots[i] = plot_activations(sample_ps[i], t_pts, 3; ymax=ymax, lw=lw, xguide="", yguide="", saveat=0.2)

# Handles special case: Stable bistability.
i = 2
plot_simulation(sample_ps[i], (0.0, t_pts[3]); u0=[0.85, 0.85, 0.85, 0.85], ymax=ymax, color=blue_scale[1], la=0.8, lw=lw, saveat=0.2)
plots[2] = plot_activation!(sample_ps[i], t_pts; idxs = [1], ymax=ymax, color=blue_scale[2], la=0.8, lw=lw, xguide="", yguide="", saveat=0.2)

# Make combined plot.
reordered_plots = [
    plots[5], plots[11], plots[10], plots[12], plots[16], 
    plots[4], plots[3], plots[2], 
    plots[15], plots[14], plots[13], 
    plots[1], plots[7], plots[6], plots[8], plots[9]]
annotated_behaviour_map = plot(reordered_plots[1:8]...,behaviour_plot, reordered_plots[9:16]...; 
     size = (2000,1200), layout = @layout [grid(1,5){0.2h}; [[a; b; c{0.28h}] d{0.595w} [e; f; g{0.28h}]]; grid(1,5){0.2h}])

# Make extra behaviour grid plots.

bg2 = BehaviourGrid([0.3,0.05,3,0.1], dataset);
bg_plot_2 = plot_behaviour_grid(bg2,start_s_slice=100,idx_axis=false,xguide="",yguide="",xticks=[],yticks=[])
bg_plot_2 = plot!(bg_plot_2, top_margin=20mm, left_margin=2mm)

bg3 = BehaviourGrid([2.0,0.2,3,0.1], dataset);
bg_plot_3 = plot_behaviour_grid(bg3,start_s_slice=100,idx_axis=false,xguide="",yguide="",xticks=[],yticks=[])
bg_plot_3 = plot!(bg_plot_3, top_margin=20mm, left_margin=2mm)

bg4 = BehaviourGrid([2.0,0.05,3,0.02], dataset);
bg_plot_4 = plot_behaviour_grid(bg4,start_s_slice=100,idx_axis=false,xguide="",yguide="",xticks=[],yticks=[])
bg_plot_4 = plot!(bg_plot_4, top_margin=20mm, left_margin=2mm)

bg5 = BehaviourGrid([2.0,0.2,3,0.01], dataset);
bg_plot_5 = plot_behaviour_grid(bg5,start_s_slice=100,idx_axis=false,xguide="",yguide="",xticks=[],yticks=[])
bg_plot_5 = plot!(bg_plot_5, top_margin=20mm, left_margin=2mm)

figure_3_full = plot(annotated_behaviour_map, bg_plot_2, bg_plot_3, bg_plot_4, bg_plot_5; size = (2000,1700), layout = @layout [a{0.77h}; grid(1,4)])

# Save plot.
make_figure(figure_3_full, "Figure_3")