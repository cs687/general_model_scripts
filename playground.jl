### Environment Setup ###

# Activate environment.
using Pkg
Pkg.activate(".")

# Add packages (first time only)
Pkg.add("Plots")
Pkg.add("MAT")

# Check which packages are in environment.
Pkg.status("MAT")

### Package Activation ###
using Plots
using MAT

### File Loading ###

### Trajectory Preparation ###

### Plotting ###