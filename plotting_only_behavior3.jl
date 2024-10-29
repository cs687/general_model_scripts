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
using Measures

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
# Define the 8x8 grid layout
# layout = @layout grid(2,2)

# Initialize an empty list to hold the plots
plots = []
tau_val= [0.1,0.2,0.5,1,2,5,10,20,50]
nu_val=[0.001,0.002,0.005,0.01,0.02,0.05,0.1]
v0_val=[0.01,0.02,0.05,0.1,0.15,0.2]
n_val=[2,3,4]
# tau_val= [0.1,1,10]
# nu_val=[0.001,0.1]
# tau_val= [0.1,50]
# nu_val=[0.001,0.1]

for v0 in eachindex(v0_val),n= eachindex(n_val)

#for n in eachindex(n_val)
    local layout = @layout grid(length(nu_val),length(tau_val))
    #local v0=1
    
    #reset plots for every iteration
    local plots=[]


    local title = plot(title = "Plot title", grid = false, showaxis = false, bottom_margin = -50Plots.px)

    for nu in eachindex(nu_val), tau in eachindex(tau_val) 
        local params = [tau_val[tau],v0_val[v0],n_val[n],nu_val[nu]]; #(τ, v0, n, η) might work: params = [0.1,0.02,3,0.1]
        local bg = BehaviourGrid(params,dataset); 
        # p=plot_behaviour_grid(bg,start_s_slice=100,idx_axis=false,xguide="D",yguide="S")
        local p=plot_behaviour_grid2(bg,start_s_slice=100,idx_axis=false)
        # tau=parse(Int, "$tau_val") 
        # nu=parse(Int, "$nu_val") 
        # Only show y-axis labels for the leftmost plots (column 1)
        if tau==1 && nu==1
            plot!(p,ylabel="nu:"*string(nu_val[nu])*"\nS",xtick=false,title="tau:"*string(tau_val[tau]))
        elseif tau==1 && nu==length(nu_val)
            plot!(p,xlabel="D\n"*"v0="*string(v0_val[v0])*" n="*string(n_val[n]),ylabel="nu:"*string(nu_val[nu])*"\nS") 
        elseif tau==1 && nu!=length(nu_val[nu])
            plot!(p,ylabel="nu:"*string(nu_val[nu])*"\nS",xtick=false)
        elseif tau!=1 && nu==length(nu_val)
            plot!(p,xlabel="D",ytick=false)
        elseif tau!=1 && nu==1
            plot!(p,ticks=false,title="tau:"*string(tau_val[tau]))
        else
            plot!(p,ticks=false)
        end 
        


        # Append the plot to the list
        push!(plots, p)



    end

    # Combine all the plots into a single figure with the 8x8 layout
    plot_now=plot(plots..., layout = layout, size=(1550, 1200),left_margin = 0mm, right_margin = 0mm)
    display(plot_now)
    savefig(joinpath("Figures","v0_"*string(v0_val[v0])*"-n_"*string(n_val[n])*".png")) 
    layout=nothing
 end