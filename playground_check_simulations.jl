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
include("../Codebase/playground_analysis.jl")

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
gr();
default(fmt = :png, framestyle=:box, grid=false, xguidefontsize=16, yguidefontsize=16, titlefontsize=18);

# Sets colour scales.
red_scale = make_red_scale(4);
bluegreen_scale = make_bluegreen_scale(4);
blue_scale= make_blue_scale(3);

##################
### Playground ###
##################


# Plot gird.
bg = BehaviourGrid([2.0,0.15,3,0.05], dataset)
plt1 = plot_behaviour_grid(bg, start_s_slice=100, idx_axis=false, xguide="", yguide="", xticks = [], yticks = [])

bg = BehaviourGrid([0.1,0.02,3,0.05], dataset)
plt2 = plot_behaviour_grid(bg, start_s_slice=100, idx_axis=false, xguide="", yguide="", xticks = [], yticks = [])

savefig(plt1, "fig1.svg")
savefig(plt2, "fig2.svg")

params = [0.1,0.02,3,0.1]
bg = BehaviourGrid(params, dataset)
plot_behaviour_grid(bg, start_s_slice=100, idx_axis=false, xguide="D", yguide="S")



plt1_1, plt1_2 = selected_and_simulated(7.0, 0.1, [0.1,0.02,3,0.05], dataset, 3; marker_color=RGB{Float64}(1.,0.2,0.2), adaptive=false, dt = 0.01)
plot(plt1_1, plt1_2; size = (1200,500))

plt1_1, plt1_2 = selected_and_simulated(1.9, 0.2, [2.0,0.15,3,0.05], dataset, 3; marker_color=RGB{Float64}(1.,0.2,0.2), adaptive=false, dt = 0.01)
plot(plt1_1, plt1_2; size = (1200,500))

plt1_1, plt1_2 = selected_and_simulated(4.2, 4.9, [0.1,0.02,3,0.1], dataset, 3; marker_color=RGB{Float64}(1.,0.2,0.2), adaptive=false, dt = 0.01)
plot(plt1_1, plt1_2; size = (1200,500))


plt2_1, plt2_2 = selected_and_simulated(2.3, 0.5, [0.1,0.02,3,0.1], dataset, 3; marker_color=RGB{Float64}(1.,0.2,0.2))
plt3_1, plt3_2 = selected_and_simulated(2.3, 0.7, [0.1,0.02,3,0.1], dataset, 3; marker_color=RGB{Float64}(1.,0.2,0.2))
plt4_1, plt4_2 = selected_and_simulated(2.3, 1.0, [0.1,0.02,3,0.1], dataset, 3; marker_color=RGB{Float64}(1.,0.2,0.2))
plt5_1, plt5_2 = selected_and_simulated(2.3, 1.3, [0.1,0.02,3,0.1], dataset, 3; marker_color=RGB{Float64}(1.,0.2,0.2))
plt6_1, plt6_2 = selected_and_simulated(2.3, 1.8, [0.1,0.02,3,0.1], dataset, 3; marker_color=RGB{Float64}(1.,0.2,0.2))
p1 = plot(plt1_1, plt2_1, plt3_1, plt4_1, plt5_1, plt6_1, plt1_2, plt2_2, plt3_2, plt4_2, plt5_2, plt6_2; layout=(2,6), size=(3000, 1200))


# Returns plots showing where a parameter set is on a map, and it's simulations.
function selected_and_simulated(S, D, params, dataset, n; start_s_slice=100, activation_lw=5, la=0.5, saveat=1.0,xticks=[], marker_color=RGB{Float64}(0.,0.,0.), ymax = 1.25, adaptive = true, dt = 0.001)    
    D_idx = argmin(abs.(dataset.D_grid .- D))
    S_idx = argmin(abs.(dataset.S_grid .- S)) - start_s_slice + 1
    bg = BehaviourGrid(params, dataset)
    plt1 = plot_behaviour_grid(bg; start_s_slice=start_s_slice, idx_axis=false, xguide="",yguide="", markers=[(S_idx, D_idx)],marker_color=marker_color)
    plt2 = plot_activations([S; D; params], (-200.0, 0.0, 1000.0), n; activation_lw=activation_lw, la=la, saveat=saveat, adaptive= adaptive, dt = dt, ymax = 1.25, xticks=xticks)
    return plt1, plt2
end

plot(plt1_1, plt1_2)

@time plt1_1, plt1_2 = selected_and_simulated(5, 50., [0.1,0.02,3,0.1], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2),xticks=[0.0, 1000.0], adaptive=false, dt = 0.01)
@time plt2_1, plt2_2 = selected_and_simulated(15, 50., [0.1,0.02,3,0.1], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2),xticks=[0.0, 1000.0], adaptive=false, dt = 0.01)
@time plt3_1, plt3_2 = selected_and_simulated(25, 50., [0.1,0.02,3,0.1], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2),xticks=[0.0, 1000.0], adaptive=false, dt = 0.01)
@time plt4_1, plt4_2 = selected_and_simulated(40, 50., [0.1,0.02,3,0.1], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2),xticks=[0.0, 1000.0], adaptive=false, dt = 0.01)
@time plt5_1, plt5_2 = selected_and_simulated(60, 50., [0.1,0.02,3,0.1], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2),xticks=[0.0, 1000.0], adaptive=false, dt = 0.01)
p1 = plot(plt1_1, plt2_1, plt3_1, plt4_1, plt5_1, plt1_2, plt2_2, plt3_2, plt4_2, plt5_2; layout=(2,5), size=(3600, 1200))

function get_mean_ss_activty(ps; n=10, l1 = 3000, l2 = 3500, saveat=1.0, adaptive = false, dt = 0.01)
    sols = cle_monte(ps, (0.0, Float64(l2)), n; saveat=saveat, adaptive = adaptive, dt = dt, u0s=fill([1.0, 1.0, 1.0, 1.0], n))
    vals = vcat([getindex.(sol.u, 1)[1000:end] for sol in sols.u]...)
    return mean(vals), std(vals)
end

function get_activation_fraction(ps; n=100, l1 = 3000, l2 = 3500, saveat=1.0, adaptive = false, dt = 0.01)
    sols = cle_monte_activation(ps, (0.0, 500.0), 250.0, n; saveat = saveat, adaptive = adaptive, dt = dt)
    return count(sol.u[end][1] > 0.5 for sol in sols)/n
end

ps = [0.1, 0.02, 3, 0.05];
n=10
l1 = 3000
l2 = 3500
saveat=1.0
adaptive = false
dt = 0.01

sols = cle_monte_activation(ps, (0.0, 500.0), 250.0, n; saveat = saveat, adaptive = adaptive, dt = dt)
frac = count(sol.u[end][1] > 0.5 for sol in sols)

@time get_activation_fraction([7.0, 0.1, 0.1, 0.02, 3, 0.05]; n = 100)

@time act_times = [get_activation_fraction([6.0, D, 0.1, 0.02, 3, 0.05]) for D in 10 .^(range(-1.0,stop=2,length=30))]
@time act_times2 = [get_activation_fraction([1.9, D, 2.0,0.15,3,0.05]) for D in 10 .^(range(-1.0,stop=2,length=30))]

act_plots = [plot_activations([S, 50.0, 0.1, 0.02, 3, 0.05], (-100.0, 0.0, 500.0), 3; ymax = 1.25) for S in 10 .^(range(0.0,stop=2,length=20))]

fig5 = plot(act_plots..., size = (2000,1200); xguide="", yguide="")
savefig(fig5, "fig5.svg")

S = (10 .^(range(0.0,stop=2,length=20)))[17]
act_plots[17] = plot_activations([S, 50.0, 0.1, 0.02, 3, 0.05], (-100.0, 0.0, 500.0), 3; ymax = 1.25)

S

plot_activations([5.0, 50.0, 0.1, 0.02, 3, 0.05], (-100.0, 0.0, 500.0), 3)
plot_activations([50.0, 50.0, 0.1, 0.02, 3, 0.05], (-100.0, 0.0, 500.0), 3)
plot_activations([50.0, 50.0, 0.1, 0.02, 3, 0.05], (-100.0, 0.0, 500.0), 3)
plot_activations([50.0, 50.0, 0.1, 0.02, 3, 0.05], (-100.0, 0.0, 500.0), 3)
plot_activations([50.0, 50.0, 0.1, 0.02, 3, 0.05], (-100.0, 0.0, 500.0), 3)




plt3 = plot(10 .^(range(-1.0,stop=2,length=20)), act_times2[1:20], xaxis=:log10; lw = 5, xticks=[], yticks = [], label = "")
savefig(plt3, "plt4.svg")
[]
mean_vals = [get_mean_ss_activty([5.0, D, 0.1, 0.02, 3, 0.05]) for D in 10 .^(range(-1.0,stop=2,length=40))]
plot()

get_mean_ss_activty([1.3, 0.5, 0.1,0.02,3,0.1])
@time get_mean_ss_activty([1.3, 0.5, 0.1,0.02,3,0.1])

D_vals = [0.1, 1.0, 3.0, 10.0, 50.0]
S_vals = 10 .^(range(-1.0,stop=2,length=40))
mean_vals = [get_mean_ss_activty([S, D, 0.1, 0.02, 3, 0.1]) for S in S_vals, D in D_vals]


isparameter()

plot(S_vals, first.(mean_vals[:,1]); ribbon = last.(mean_vals[:,1]), xaxis=:log10, lw = 4, label="D = $(D_vals[1])", xlimit = (S_vals[1], S_vals[end]), legendfontsize=12)
plot!(S_vals, first.(mean_vals[:,2]); ribbon = last.(mean_vals[:,2]), xaxis=:log10, lw = 4, label="D = $(D_vals[2])", xlimit = (S_vals[1], S_vals[end]), legendfontsize=12)
plot!(S_vals, first.(mean_vals[:,3]); ribbon = last.(mean_vals[:,3]), xaxis=:log10, lw = 4, label="D = $(D_vals[3])", xlimit = (S_vals[1], S_vals[end]), legendfontsize=12)
plot!(S_vals, first.(mean_vals[:,4]); ribbon = last.(mean_vals[:,4]), xaxis=:log10, lw = 4, label="D = $(D_vals[4])", xlimit = (S_vals[1], S_vals[end]), legendfontsize=12)
plot!(S_vals, first.(mean_vals[:,5]); ribbon = last.(mean_vals[:,5]), xaxis=:log10, lw = 4, label="D = $(D_vals[5])", xlimit = (S_vals[1], S_vals[end]), legendfontsize=12)
p1 = plot!(xguide = "S value", yguide = "Steady state \nsigma activity (mean +-std)")

plot(S_vals, first.(mean_vals[:,1]); ribbon = last.(mean_vals[:,1]), xaxis=:log10, lw = 4, label="D = $(D_vals[1])", xlimit = (S_vals[1], S_vals[end]), legendfontsize=12)

plt3 = plot(10 .^(range(-1.0,stop=2,length=40)), first.(mean_vals), xaxis=:log10; lw = 5, xticks=[], yticks = [], label = "")
savefig(plt3, "plt3.svg")

S_vals[12]

params = [0.1,0.02,3,0.1]
bg = BehaviourGrid(params, dataset)
plot_behaviour_grid(bg, start_s_slice=100, idx_axis=false, xguide="D", yguide="S")

plot!([2, 2], [1, 201], lw = 4, color = 1, label="D = $(D_vals[1])")
plot!([findfirst(D_grid .> D_vals[2]), findfirst(D_grid .> D_vals[2])], [1, 201], lw = 4, color = 2, label="D = $(D_vals[2])")
plot!([findfirst(D_grid .> D_vals[3]), findfirst(D_grid .> D_vals[3])], [1, 201], lw = 4, color = 3, label="D = $(D_vals[3])")
plot!([findfirst(D_grid .> D_vals[4]), findfirst(D_grid .> D_vals[4])], [1, 201], lw = 4, color = 4, label="D = $(D_vals[4])")
p2 = plot!([findfirst(D_grid .> D_vals[5]), findfirst(D_grid .> D_vals[5])], [1, 201], lw = 4, color = 5, label="D = $(D_vals[5])")

plot(p1 ,p2, size = (1200,400), bottom_margin = 10mm, left_margin = 15mm)


D_grid

ps = [1.3, 0.5, 0.1,0.02,3,0.1]
cle_monte(p::Vector{Float64}, tspan::Tuple{Float64,Float64}, n::Int64; model=model3, u0_fp=1, u0s=get_cle_u0s(u0_fp, model, p, n), callbacks=[], kwargs...)

ps = [1.3, 0.5, 0.1,0.02,3,0.1]
sols = cle_monte(ps, (0.0, 2000.0), 5, saveat=1.0, adaptive = false, dt = 0.01, u0s=fill([1.0, 1.0, 1.0, 1.0], 5))
vcat([getindex.(sol.u, 1)[1000:end] for sol in sols.u]...)


plot(sols)

plot(plt4_1, plt4_2)

plt11 = plot_activation([7.5, 10.0, 0.1,0.02,3,0.05], (-200.0, 0.0, 1000.0); ymax = 1.25, idxs = [1], title = "S = 7.5, D = 10.0")
plt21 = plot_activation([15.0, 10.0, 0.1,0.02,3,0.05], (-200.0, 0.0, 1000.0); ymax = 1.25, idxs = [1], title = "S = 15.0, D = 10.0")
plt12 = plot_activation([7.5, 20.0, 0.1,0.02,3,0.05], (-200.0, 0.0, 1000.0); ymax = 1.25, idxs = [1], title = "S = 7.5, D = 20.0")
plt22 = plot_activation([15.0, 20.0, 0.1,0.02,3,0.05], (-200.0, 0.0, 1000.0); ymax = 1.25, idxs = [1], title = "S = 15.0, D = 20.0")

plot(plt21, plt22, plt11, plt12, layout = (2,2), size = (1200,700), bottom_margin = 10mm)


plot_activation([5, 0.01, 0.1,0.02,3,0.1], (-200.0, 0.0, 1000.0); ymax = 1.25, idxs = [1])
plot_activation([10, 0.01, 0.1,0.02,3,0.1], (-200.0, 0.0, 1000.0); ymax = 1.25, idxs = [1])


plot_behaviour_grid(bg, start_s_slice=100, idx_axis=false, xguide="D", yguide="S")
plot_activation([10, 20.01, 0.1,0.02,3,0.1], (-200.0, 0.0, 1000.0); ymax = 1.25, idxs = [1])
plot_activation([75, 20.01, 0.1,0.02,3,0.1], (-200.0, 0.0, 1000.0); ymax = 1.25, idxs = [1])

plot_activation([10, 25.01, 0.1,0.02,3,0.1], (-200.0, 0.0, 1000.0); ymax = 1.25, idxs = [1])
plot_activation([15, 25.01, 0.1,0.02,3,0.1], (-200.0, 0.0, 1000.0); ymax = 1.25, idxs = [1])

plt1_1, plt1_2 = selected_and_simulated(2.3, 0.2, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
plot(plt1_1, plt1_2; layout = (2,1), size=(1200,600))
plot(plt1_2)


### Old Stuff ###

plt1_1, plt1_2 = selected_and_simulated(2.3, 0.2, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
plt2_1, plt2_2 = selected_and_simulated(2.3, 0.5, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
plt3_1, plt3_2 = selected_and_simulated(2.3, 0.7, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
plt4_1, plt4_2 = selected_and_simulated(2.3, 1.0, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
plt5_1, plt5_2 = selected_and_simulated(2.3, 1.3, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
plt6_1, plt6_2 = selected_and_simulated(2.3, 1.8, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
p1 = plot(plt1_1, plt2_1, plt3_1, plt4_1, plt5_1, plt6_1, plt1_2, plt2_2, plt3_2, plt4_2, plt5_2, plt6_2; layout=(2,6), size=(3000, 1200))

plt1_1, plt1_2 = selected_and_simulated(2.8, 0.2, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
plt2_1, plt2_2 = selected_and_simulated(2.8, 0.5, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
plt3_1, plt3_2 = selected_and_simulated(2.8, 0.7, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
plt4_1, plt4_2 = selected_and_simulated(2.8, 1.0, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
plt5_1, plt5_2 = selected_and_simulated(2.8, 1.3, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
plt6_1, plt6_2 = selected_and_simulated(2.8, 1.8, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
plt7_1, plt7_2 = selected_and_simulated(2.8, 3, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
plt8_1, plt8_2 = selected_and_simulated(2.8, 6, [2.0,0.1,3,0.05], dataset, 5; marker_color=RGB{Float64}(1.,0.2,0.2))
p2 = plot(plt1_1, plt2_1, plt3_1, plt4_1, plt5_1, plt6_1, plt7_1, plt8_1, plt1_2, plt2_2, plt3_2, plt4_2, plt5_2, plt6_2, plt7_2, plt8_2; layout=(2,8), size=(4000, 1200))

plt1_1, plt1_2 = selected_and_simulated(2.3, 0.12, [10.0,0.05,3,0.1], dataset, 2; marker_color=RGB{Float64}(1.,0.2,0.2))
plt2_1, plt2_2 = selected_and_simulated(2.3, 0.24, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
plt3_1, plt3_2 = selected_and_simulated(2.3, 0.48, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
p1 = plot(plt1_1, plt2_1, plt3_1, plt1_2, plt2_2, plt3_2; layout=(2,3), size=(1800, 1200))

plt1_1, plt1_2 = selected_and_simulated(2.45, 0.12, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
plt2_1, plt2_2 = selected_and_simulated(2.45, 0.24, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
plt3_1, plt3_2 = selected_and_simulated(2.45, 0.48, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
p2 = plot(plt1_1, plt2_1, plt3_1, plt1_2, plt2_2, plt3_2; layout=(2,3), size=(1800, 1200))

plt1_1, plt1_2 = selected_and_simulated(3.0, 0.12, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
plt2_1, plt2_2 = selected_and_simulated(3.0, 0.24, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
plt3_1, plt3_2 = selected_and_simulated(3.0, 0.48, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
p3 = plot(plt1_1, plt2_1, plt3_1, plt1_2, plt2_2, plt3_2; layout=(2,3), size=(1800, 1200))

plt1_1, plt1_2 = selected_and_simulated(3.65, 0.12, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
plt2_1, plt2_2 = selected_and_simulated(3.65, 0.24, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
plt3_1, plt3_2 = selected_and_simulated(3.65, 0.48, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
p4 = plot(plt1_1, plt2_1, plt3_1, plt1_2, plt2_2, plt3_2; layout=(2,3), size=(1800, 1200))

plt1_1, plt1_2 = selected_and_simulated(3.95, 0.12, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
plt2_1, plt2_2 = selected_and_simulated(3.95, 0.24, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
plt3_1, plt3_2 = selected_and_simulated(3.95, 0.48, [10.0,0.05,3,0.1], dataset, 50; marker_color=RGB{Float64}(1.,0.2,0.2))
p5 = plot(plt1_1, plt2_1, plt3_1, plt1_2, plt2_2, plt3_2; layout=(2,3), size=(1800, 1200))

plot(p4)


