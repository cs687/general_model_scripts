# Define custom ticks for the y-axis (logarithmic scale)

x=1:5
y = 10 .^ x  # 10 raised to the power of x, i.e., 10^1, 10^2, ..., 10^5

yticks = [10^i for i in 1:5]  # Ticks at 10^1, 10^2, 10^3, etc.
# Create the plot with custom y-ticks
plot(x, y, yscale=:log10, yticks=(yticks, string.(yticks)), xlabel="x", ylabel="log scale (y)")

# Show the plot
#display()
