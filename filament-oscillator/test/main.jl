using DifferentialEquations
using ODEInterfaceDiffEq

include("../src/models/cilia.jl")
include("../src/models/solver.jl")
include("../src/models/physical_params.jl")


# Initial phase
Ψ₀ = repeat([π + 0.3, 0.0], outer=[M])

# Simulation parameters
num_periods = 5
alg = Trapezoid(autodiff=false)
# alg = Euler()
# Number of timesteps. Set to 0 if using adaptive timestepping
# num_steps = 1001
num_steps = 0

# Instance objects
h = L/(N - 1)
params = CiliaParameters(M, N, φ, d, a, κ_b)
x₀ = [[(j - 1)*d, 0.0, 0.0] for j=1:M]
s = [(i)*h for i=1:N]
A = [cos(φ) -sin(φ) 0.0; sin(φ) cos(φ) 0.0; 0.0 0.0 1.0]
system = CiliaSystem(
    params,
    x₀,
    s,
    A,
    Ψ₀
)
fluid = FluidParameters(μ)

# Run the system
solution = run_system(system, fluid, num_periods*T, alg, num_steps)
