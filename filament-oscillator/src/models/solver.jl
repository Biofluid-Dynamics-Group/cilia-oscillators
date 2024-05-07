using LinearSolve
using DifferentialEquations

include("cilia.jl")


"""
    Ψ_dot(system::CiliaSystem, fluid::FluidParameters)

Returns the phase velocity at a given system phase. The phase velocity is determined by
solving the saddle-point system that serves as the equation of motion for the system of
filaments.
"""
function Ψ_dot(system::CiliaSystem, fluid::FluidParameters)
    Q = Q_ref(system, fluid)
    κ_matrix = κₕ(system)
    matrix = [
        Πₕ(system, fluid) -κ_matrix;
        -κ_matrix' zeros(2system.params.M, 2system.params.M)
    ]
    elastic_relaxation = repeat([0.0, system.κ_b*system.Ψ[2j]]; outer=[system.params.M])
    force = Q + elastic_relaxation
    rhs = vcat(zeros(system.M*3system.N), -force)
    problem = LinearProblem(matrix, rhs)
    solution = solve(problem, KrylovJL_GMRES())
    dΨ_dt = solution[end-(2*system.M - 1):end]
    return dΨ_dt
end

"""
    filament_oscillators!(dΨ_dt::Vector, Ψ::Vector, p::NamedTuple, t::Real)

Function to construct the ODEProblem object. Updates the phase velocities `dΨ_dt` given
the current phase `Ψ`. `p` holds the system. Since the system is autonomous, the variable
`t` is not used.
"""
function filament_oscillators!(dΨ_dt::Vector, Ψ::Vector, p::NamedTuple, t::Real)
    # The solver tells the system its new phase
    p.system.Ψ .= Ψ
    # The system tells the solver the new phase velocities given the current phase
    dΨ_dt .= Ψ_dot(system, p.fluid)
    return nothing
end

"""
    run_system(
        system::CiliaSystem, fluid::FluidParameters, time::Real,
        alg::Union{DEAlgorithm,Nothing}
    )

Given a system of cilia in its current state, advance the system for `time` milliseconds
using the algorithm `arg` for time integration.
"""
function run_system(system::CiliaSystem, fluid::FluidParameters, time::Real, alg)
    t_span = (0.0, time)
    Ψ₀ = copy(system.Ψ)
    p = (system=system, fluid=fluid)
    problem = ODEProblem(filament_oscillators!, Ψ₀, t_span, p)
    solution = solve(problem, alg)
    return solution
end
