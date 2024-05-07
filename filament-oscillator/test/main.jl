include("../models/multiple_cilia.jl")


# Physical parameters
a = 7e-2  # μm
μ = 1.0  # Pa s
κ_b = 1.0  # pN μm^2
T = 60  # ms

# Simulation parameters
num_steps = 500
num_periods = 2
alg = radau()

# Run the system
p = (a=a, μ=μ, κ_b=κ_b)
solution = run_system(p, num_periods*T, num_steps, alg)

# Plot the system trajectory
plot_system_trajectory(solution)
