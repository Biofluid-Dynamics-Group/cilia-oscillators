include("../src/models/cilia.jl")
include("../src/visualisation/plotters.jl")


# Physical parameters
a = 7e-2  # μm, sphere radius
μ = 1.0  # Pa s, fluid viscosity
κ_b = 1.0  # pN μm^2, bending stiffness

# Geometric setup
M = 1  # Number of cilia
d = 40.0  # μm, distance between cilia
φ = 12.5*π/180.0  # rad, angle of cilia beat plane

# Discretisation parameters
N = 5  # Number of discretisation spheres per cilium

# Initial phase
Ψ₀ = repeat([2π, 0.0], outer=[M])

# Simulation parameters
num_periods = 1
alg = radau()

# Instance objects
h = L/(N - 1)
system = CiliaSystem(
    params = CiliaParameters(M, N, φ, d, a, κ_b),
    x₀ = [[(j - 1)*d, 0.0, 0.0] for j=1:M],
    s = [(i)*h for i=1:N]
    A = [cos(φ) -sin(φ) 0.0; sin(φ) cos(φ) 0.0; 0.0 0.0 1.0],
    Ψ = Ψ₀
)
fluid = FluidParameters(μ)

# Run the system
solution = run_system(system, fluid, num_periods*T, alg)
