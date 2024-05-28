include("src/utils/integral.jl")
using .Integral
using DifferentialEquations
using ODEInterfaceDiffEq

include("src/models/solver.jl")
include("src/models/physical_params.jl")


# Initial phase
Ψ₀ = repeat([π + 0.3, 0.0], outer=[M])

# Simulation parameters
num_periods = 1
# alg = Trapezoid(autodiff=false)
alg = ABM54()
# Number of timesteps. Set to 0 if using adaptive timestepping
num_steps = 1000
# num_steps = 0

# Define background flow
function u(t::Real)
    return [0., 0., 0.]
end

# Integration rule
rule = GaussLegendre(15)

# Instance objects
beat_params = BeatParameters(L, T, θ_0, f_eff, f_ψ, f_w, orientation, rule)
h = L/(N - 1)
sim_params = SimulationParameters(M, N, φ, d, a, κ_b)
x₀ = [[(j - 1)*d, 0.0, 0.0] for j=1:M]
s = [(i)*h for i=1:N]
A = [cos(φ) -sin(φ) 0.0; sin(φ) cos(φ) 0.0; 0.0 0.0 1.0]
system = CiliaSystem(
    beat_params,
    sim_params,
    x₀,
    s,
    A,
    Ψ₀,
    u
)
fluid = FluidParameters(μ)

# Run the system
solution = run_system(system, fluid, num_periods*T, alg, num_steps)
