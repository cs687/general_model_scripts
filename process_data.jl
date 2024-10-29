### Get Behaviour Classifications ###
# From a behaviour mapping, finds and saves all behaviour_classifications.

### Fetch Functions ###
include("../Codebase/file_management.jl")
include("../Codebase/behaviour_grid_loading.jl")
include("../Codebase/basic_analysis.jl")
include("../Codebase/models.jl")
include("../Codebase/model_simulation.jl")
include("../Codebase/utility.jl")

### Declares Behaviour Map ###
S_grid = 10 .^(range(-1,stop=2,length=300))
D_grid = 10 .^(range(-1,stop=2,length=300))
τ_grid = [0.1,0.15,0.20,0.30,0.50,0.75,1.0,1.5,2.0,3.0,5.0,7.50,10.0,15.0,20.0,30.0,50.0,75.0,100.0]
v0_grid = [0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.20]
n_grid = [2.0,3.0,4.0]
η_grid = [0.001,0.002,0.005,0.01,0.02,0.05,0.1]

target_folder = "BehaviourMap"
dataset = DataSet(target_folder,S_grid,D_grid,τ_grid,v0_grid,n_grid,η_grid);

### Behaviour Map Parameters ###
# Saves the information of what behaviour each parameter set was assigned to (in a compressed form).
let
    behaviours = fill(:unfilled, length.([S_grid, D_grid, τ_grid, v0_grid, n_grid, η_grid])...)
    for (τ_idx,τ) in enumerate(τ_grid), (v0_idx,v0) in enumerate(v0_grid), (n_idx,n) in enumerate(n_grid), (η_idx,η) in enumerate(η_grid)
        bg = BehaviourGrid([τ, v0, n, η], dataset)
        for (S_idx,S) in enumerate(S_grid), (D_idx,D) in enumerate(D_grid)
            behaviours[S_idx ,D_idx ,τ_idx, v0_idx, n_idx, η_idx] = bg.behaviours[S_idx, D_idx]
        end
    end
    
    serialize("Data/Simulated_data/$(dataset.dataset_tag)/all_behaviour_classifications.jls", behaviours)
end


### Simulation Examples ###

# Simulated examples for specific parameter values.
let 
    n_transitions = 50
    n_samples = 100
    n_repeat = 100
    sim_leng = 2000.0

    simulations = Dict{Vector{Float64}, Vector}()
    D = 10.0
    D_idx = findfirst(D_grid .> D)
    τ = 1.0; v0 = 0.05; n = 3; η = 0.05;

    bg = BehaviourGrid([τ, v0, n, η], dataset)

    sim_results = Vector{Vector{Vector{Float64}}}()
    for S in logrange(S_grid[101], S_grid[end], n_samples)
        sols = cle_monte([S, D, τ, v0, n, η], sim_leng, 2*n_repeat; saveat=0.2, u0=sys_roots(model3,[S, D, τ, v0, n, η])[end]).u
        vector_sols = [sol[1,:] for sol in sols]
        filtered_sols = filter(sol -> length(sol) == 10001, vector_sols)
        (length(filtered_sols) < 100) && (@warn "Few simulations ($(length(filtered_sols))) for S=$S.")
        push!(sim_results, filtered_sols[1:100])
    end

    serialize("Data/Simulated_data/$(dataset.dataset_tag)/transition_simulation_known_example.jls", sim_results)
end

sim_results = deserialize("Data/Simulated_data/$(dataset.dataset_tag)/transition_simulation_known_example.jls")

# Creates and saves example trajectories for the combined simulation/experiment classifier to work on.
let 
    n_transitions = 50
    n_samples = 100
    n_repeat = 100
    sim_leng = 2000.0

    simulations = Dict{Vector{Float64}, Vector}()

    @time while counts < 100
        D_idx = rand(1:length(D_grid)); D = D_grid[D_idx];
        τ = rand(τ_grid); v0 = rand(v0_grid); n = rand(n_grid); η = rand(η_grid[1:end-1]);
        bg = BehaviourGrid([τ, v0, n, η], dataset)
        (count(bg.behaviours[:,D_idx] .== :no_activation) < 10) && continue
        (count(bg.behaviours[:,D_idx] .== :stochastic_pulsing) < 10) && continue
        (count(bg.behaviours[:,D_idx] .== :oscillation) < 10) && continue
        (count(bg.behaviours[:,D_idx] .== :homogeneous_activation) < 10) && continue
        println(counts)

        sim_results = [[sol[1,:] for sol in cle_monte([S, D, τ, v0, n, η], sim_leng, n_repeat; saveat=1.0, u0=sys_roots(model3,[S, D, τ, v0, n, η])[end]).u] for S in logrange(S_grid[101], S_grid[end], n_samples)]
        all(length.(vcat(sim_results...)) .== 2001) || continue
        simulations[[D, τ, v0, n, η]] = sim_results
        counts += 1
    end

    serialize("Data/Simulated_data/$(dataset.dataset_tag)/transition_simulation_examples.jls", simulations)
end