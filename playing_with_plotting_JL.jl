using Plots

# Define the 8x8 grid layout
layout = @layout grid(8, 7)

# Initialize an empty list to hold the plots
plots = []

# Generate the 8x8 grid of subplots using a for loop
for i in 1:8, j in 1:7
    # Create a simple line plot with random data
    
    p = plot(1:10, rand(10), xlims=[0,10], ylims=[0,1],legend=false)
    # Only show y-axis labels for the leftmost plots (column 1)
    if i==1 && j==1
        plot!(p,ylabel="MY",xtick=false,title="pups2")
    elseif j==1 && i!=8
        plot!(p,ylabel="MY",xtick=false)
    elseif j!=1 && i==8
        plot!(p,xlabel="time",ytick=false)
    elseif j==1 && i==8
        plot!(p,xlabel="time",ylabel="MY")
    elseif i==1 && j!=1
        plot!(p,ticks=false,title="pups")
    else
        plot!(p,ticks=false)
    end


    # Append the plot to the list
    push!(plots, p)
end

# Combine all the plots into a single figure with the 8x8 layout
plot(plots..., layout = layout, size=(1200, 1200))
