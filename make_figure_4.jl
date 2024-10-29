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
include("../Codebase/Functions_specific/trajectory_classification.jl")

### Behaviour Map Parameters ###

# All the parameter values of the behaviour map grid.
# 10^-1 to 10^2
S_grid = 10 .^(range(-1, stop=2, length=300))
D_grid = 10 .^(range(-1, stop=2, length=300))
τ_grid = [0.1,0.15,0.20,0.30,0.50,0.75,1.0,1.5,2.0,3.0,5.0,7.50,10.0,15.0,20.0,30.0,50.0,75.0,100.0]
v0_grid = [0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.20]
n_grid = [2.0,3.0,4.0]
η_grid = [0.001,0.002,0.005,0.01,0.02,0.05,0.1]

target_folder = "BehaviourMap"
dataset = DataSet(target_folder,S_grid,D_grid,τ_grid,v0_grid,n_grid,η_grid);


### Plot Settings ###

# Default values.
gr()
default(fmt = :png, framestyle=:box, grid=false, xguidefontsize=36, yguidefontsize=36, titlefontsize=18);

# Sets colour scales.
red_scale = make_red_scale(4);
bluegreen_scale = make_bluegreen_scale(4);
blue_scale= make_blue_scale(3);


####################
### Make Figure ####
####################

# Make behaviour grid plot.
# `params` is the parameter set [τ, v0, n, η] (full map is [S, D, τ, v0, n, η])
params = [1.0, 0.05, 3, 0.05]; # Sets parameter values.
bg = BehaviourGrid(params, dataset); # Loads the behaviour map for this parameter set. 
plot_behaviour_grid(bg, start_s_slice=100, start_d_slice=100, idx_axis=false, xguide="D", yguide="S") # Plot behaviour map.
bg_plot = plot!(size=(1800,1400))

# Makes behaviour classification plots.

# Extracts trajectories.
exp_trajss_all = load_all_trajectories()
exp_trajss = [[filter(e -> length(e) == 572, etjs) for etjs in exp_trajs] for exp_trajs in exp_trajss_all]
sim_trajs = deserialize("Data/Simulated_data/$(dataset.dataset_tag)/transition_simulation_known_example.jls")

# Plot the classification sketch (experiments).
exp_classification_plot = let
    class_fracs = [get_trajs_class_fractions(classify_trajectories(exp_trajss[i])) for i in 1:5]
    plot()
    colors = [:grey, :gold, :red3, :blue, :black]
    iptg_levels = 2:10
    for behaviour in 1:5
        means = []
        stds = []
        fracs = getindex.(class_fracs, behaviour)
        for iptg_level in iptg_levels
            rel_fracs = filter(!isnan, getindex.(fracs, iptg_level + 1))
            push!(means, mean(rel_fracs))
            push!(stds, std(rel_fracs))
        end
        plot!(iptg_levels, means; ribbon=stds, color=colors[behaviour], lw=20, la=0.8, fillalpha=0.3, xlimit=(2.0,10.0), ylimit=(-0.01, 1.01))
    end
    plot!(xguide="IPTG concentration (µM)", yguide="Behaviour fractions", bottom_margin=20mm, left_margin=30mm, right_margin=5mm)
end

# Plot the classification sketch (simulations).
sim_classification_plot = plot_classification(sim_trajs; lw=20, xguide="S (au)", yguide="Behaviour fractions")

# Function for making plot of 4 trajectories.
function plot_trajectory_examples(trajs; xguide="Time (au)", yguide="σ (au)", kwargs...)
    (length(trajs) != 4) && error("Not the right number of trajectories provided.")
    p1 = plot(trajs[1]; kwargs...)
    p2 = plot(trajs[2]; kwargs...)
    p3 = plot(trajs[3]; kwargs...)
    p4 = plot(trajs[4]; xguide=xguide, kwargs...)
    plot(p1,p2,p3,p4, layout=(4,1), xlimit=(0.0, length(trajs[1])))
    plot!(yguide=yguide)
end

# Plot simulated example trajectories.
get_sim_trajs(D_idx) = [traj[1700:end] for traj in rand(sim_trajs[D_idx], 4)]
plts_na = plot_trajectory_examples(get_sim_trajs(20); ylimit=(-0.05, 1.25), lw=8)
plts_sp = plot_trajectory_examples(get_sim_trajs(30); ylimit=(-0.05, 1.25), lw=8)
plts_o = plot_trajectory_examples(get_sim_trajs(45); ylimit=(-0.05, 1.25), lw=8)
plts_ha = plot_trajectory_examples(get_sim_trajs(75); ylimit=(-0.05, 1.25), lw=8)
sim_plot_examples = plot(plts_na, plts_sp, plts_o, plts_ha, layout=(1,4),size=(6000,1500))


# Plot experimental example trajectories.
exp = exp_trajss[2]

vals = vcat(vcat(exp...)...)
dens = kde(vals)
plot(dens.x, dens.density)


pk_idxs, pk_vals = findmaxima(dens.density)
idx_1, idx_2 = pk_idxs[sortperm(-pk_vals)][1:2]


pk_idxs, pk_vals = findmaxima(dens.density)
idx_1, idx_2 = pk_idxs[sortperm(-pk_vals)][1:2]
return dens.x[idx_1], dens.x[idx_2]

na_idxs = [1,2,3,4]
na_trajs = [exp[3][idx] for idx in na_idxs]
plts_na = plot_trajectory_examples(na_trajs;  ylimit=(-0.05, 5000), lw=8, xguide="Time (10 minutes)", yguide="σB (au)", color=:green)

sp_idxs = [10,20,30,40]
sp_trajs = [exp[5][idx] for idx in sp_idxs]
plts_sp = plot_trajectory_examples(sp_trajs;  ylimit=(-0.05, 5000), lw=8, xguide="Time (10 minutes)", yguide="σB (au)", color=:green)

o_idxs = [8,15,17,25]
o_trajs = [exp[6][idx] for idx in o_idxs]
plts_o = plot_trajectory_examples(o_trajs;  ylimit=(-0.05, 5000), lw=8, xguide="Time (10 minutes)", yguide="σB (au)", color=:green)

ha_idxs = [1,2,3,4]
ha_trajs = [exp[8][idx] for idx in ha_idxs]
plts_ha = plot_trajectory_examples(ha_trajs;  ylimit=(-0.05, 5000), lw=8, xguide="Time (10 minutes)", yguide="σB (au)", color=:green)

exp_plot_examples = plot(plts_na, plts_sp, plts_o, plts_ha, layout=(1,4),size=(6000,1500))

# Combined plot.
gr()
plt_center = plot(sim_classification_plot, bg_plot, exp_classification_plot, layout=(1,4), size=(6000,1500), left_margin=60mm, legend=:none)
figure_4 = plot(sim_plot_examples, plt_center, exp_plot_examples, size=(8000,5000), layout=@layout [a{0.28h}; b{0.44h}; c{0.28h}])


include("../Codebase/Functions_specific/trajectory_classification.jl")

# Playground stuff.
get_info = (t,s)->""
plot_classification_examples(exp_trajss[1]; info_f=get_info)
plot_classification_examples(exp_trajss[2]; info_f=get_info)
plot_classification_examples(exp_trajss[3]; info_f=get_info)
plot_classification_examples(exp_trajss[4]; info_f=get_info)
plot_classification_examples(exp_trajss[5]; info_f=get_info)

sum(length.(exp_trajss[4]))