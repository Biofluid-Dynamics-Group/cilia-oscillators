using LinearSolve
using DifferentialEquations
using DataFrames

include("cilia.jl")


"""
    Ψ_dot(system::CiliaSystem, fluid::FluidParameters)

Returns the phase velocity at a given system phase. The phase velocity is determined by
solving the saddle-point system that serves as the equation of motion for the system of
filaments.
"""
function Ψ_dot(system::CiliaSystem, fluid::FluidParameters, t::Real)
    Q = Q_ref(system, fluid)
    κ_matrix = -κₕ(system)
    Π_matrix = Πₕ(system, fluid)
    matrix = [
        Π_matrix κ_matrix;
        κ_matrix' zeros(2system.sim_params.M, 2system.sim_params.M)
    ]
    elastic_relaxation = cat(
        [[0.0, -system.sim_params.κ_b*system.Ψ[2j]] for j=1:system.sim_params.M]..., dims=1
    )
    background_flow = repeat(system.u(t), outer=[system.sim_params.M*system.sim_params.N])
    force = Q + elastic_relaxation
    rhs = vcat(background_flow, -force)
    problem = LinearProblem(matrix, rhs)
    solution = solve(problem, KrylovJL_GMRES())
    dΨ_dt = solution[end-(2*system.sim_params.M - 1):end]
    A_matrix = κ_matrix'*inv(Π_matrix)*κ_matrix
    return dΨ_dt, force, eigen(A_matrix).values, det(A_matrix)
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
    derivative, forcing, A_matrix_eigenvals, det = Ψ_dot(system, p.fluid, t)
    push!(p.derivative, derivative)
    push!(p.forcings, forcing)
    push!(p.eigenvalues, A_matrix_eigenvals)
    push!(p.determinant, det)
    dΨ_dt .= derivative
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
function run_system(
    system::CiliaSystem, fluid::FluidParameters, time::Real, alg, num_steps::Int=0
)
    t_span = (0.0, time)
    Ψ₀ = copy(system.Ψ)
    df = DataFrame()
    t = typeof(time)
    p = (
        system=system, fluid=fluid, derivative=Vector{t}[],
        eigenvalues=Vector{t}[], determinant=t[], forcings=Vector{t}[]
    )
    problem = ODEProblem(filament_oscillators!, Ψ₀, t_span, p)
    if num_steps != 0
        solution = solve(
            problem, alg=alg, adaptive=false, dt=time/num_steps, progress=true
        )
    else
        solution = solve(problem, alg=alg, progress=true)
    end
    df.time = solution.t
    Ψ = stack(solution.u, dims=1)
    for j=1:2*system.sim_params.M
        df[!, "Ψ_$j"] = Ψ[:, j]
    end
    df.dΨ = p.derivative
    df.forcings = p.forcings
    df.eigenvalues = p.eigenvalues
    df.determinant = p.determinant
    return df, solution
end
