using LinearAlgebra
using LinearSolve
using DifferentialEquations
using Plots


include("TangentAngleBeat.jl")
include("../stokes/GreensFunctions.jl")

μ = 1.0
a = 7e-2

N = 15
h = L/(N - 1)
s = [(i)*h for i=1:N]

"""
    K(s::Real, ψ::Vector)

Returns the matrix that maps the phase velocity to the filament velocity at arclength `s`
and phase `ψ`.
"""
function K(s::Real, ψ::Vector)
    return [∂ξ_∂ψ_1(s, ψ) ∂ξ_∂ψ_2(s, ψ)]
end

"""
    Kₕ(ψ::Vector)

Returns the matrix that maps the phase velocity to the discretised filament velocity at
phase `ψ`.
"""
function Kₕ(ψ::Vector)
    return cat([K(s[i], ψ) for i = 1:N]..., dims=1)
end

"""
    Mₕ(ψ::Vector, μ::Real, a::Real)

Returns the regularised Stokes mobility matrix mapping forces in the discretised filament
to velocities at phase `ψ`, in a fluid with viscosity `μ`, where the filament was discretised using spheres of radius `a`.
"""
function Mₕ(ψ::Vector, μ::Real, a::Real)
    mobility = zeros(3*N, 3*N)
    for i = 1:N
        for j = 1:N
            α = 1 + 3*(i - 1)
            β = 1 + 3*(j - 1)
            mobility[α:(α + 2), β:(β + 2)] .= rotne_prager_blake_tensor(
                    ξ(s[i], ψ), ξ(s[j], ψ), μ, a
                )
        end
    end
    return mobility
end

"""
    q_ref(ψ::Vector, μ::Real, a::Real, t::Real)

Returns the generalised force vector that induces a reference beat in a single filament.
"""
function q_ref(ψ::Vector, μ::Real, a::Real)
    return Kₕ(ψ)'*(Mₕ(ψ, μ, a)\Kₕ(ψ)*[2π/T, 0.0])
end

"""
    ψ_dot(ψ::Vector, μ::Real, a::Real, κ::Real, t::Real)

Returns the phase velocity at a given phase `ψ`, in a fluid with viscosity `μ`, where the
filament was discretised using spheres of radius `a`, and the filament has generalised
bending stiffness `κ`. The phase velocity is determined by solving the saddle-point system
that serves as the equation of motion for the filament.
"""
function ψ_dot(ψ::Vector, μ::Real, a::Real, κ_b::Real)
    q = q_ref(ψ, μ, a) .- [0.0, -κ_b*ψ[2]]
    K_matrix = Kₕ(ψ)
    matrix = [Mₕ(ψ, μ, a) -K_matrix; -K_matrix' zeros(2, 2)]
    rhs = vcat(zeros(3*N), -q)
    problem = LinearProblem(matrix, rhs)
    solution = solve(problem, SimpleGMRES())
    dψ_dt = solution[end-1:end]
    return dψ_dt
end

"""
    filament_oscillator!(dψ_dt::Vector, ψ::Vector, p::NamedTuple, t::Real)

Updates the phase velocity `dψ_dt` at phase `ψ` and time `t` with parameters `p` for
use with DifferentialEquations.jl solvers.
"""
function filament_oscillator!(dψ_dt::Vector, ψ::Vector, p::NamedTuple, t::Real)
    μ = p.μ
    a = p.a
    κ_b = p.κ_b
    rhs = ψ_dot(ψ, μ, a, κ_b)
    dψ_dt[1] = rhs[1]
    dψ_dt[2] = rhs[2]
    return nothing
end

"""
    run_filament(p::NamedTuple, final_time::Real, num_steps::Int)

Given parameters `p`, solves the filament oscillator model from t=0 to `final_time` with
`num_steps` timesteps using 4th order Adams-Bashforth.
"""
function run_filament(p::NamedTuple, final_time::Real, num_steps::Int)
    t_span = (0.0, final_time)
    problem = ODEProblem(filament_oscillator!, [2π, 0.0], t_span, p)
    solution = solve(problem, AB4(), dt=final_time/num_steps)
    return solution
end

"""
    plot_solution(solution::ODESolution)

Given an ODESolution `solution` of the filament oscillator, plots the filament's movement.
"""
function plot_solution(solution::ODESolution)
    ψ_array = stack(solution.u, dims=1)
    p = plot(
        xlim=(-L*1.1, L*1.1), ylim=(0.0, L*1.1), title="Filament movement",
        xaxis="x [μm]", yaxis="z [μm]", legend=false
    )
    color_scheme = palette(:twilight, size(ψ_array)[1], rev=true)
    for (i, ψ) in enumerate(ψ_array[:, 1])
        positions = zeros(N+1, 3)
        for j=1:N
            positions[j+1, :] += ξ(s[j], ψ_array[i, :])
        end
        plot!(positions[:, 1], positions[:, 3], color=color_scheme[i])
    end
    display(p)
    return nothing
end

"""
    plot_trajectory(solution::ODESolution)

Given an ODESolution `solution` of the filament oscillator, plots the phase evolution.
"""
function plot_trajectory(solution::ODESolution)
    ψ_array = stack(solution.u, dims=1)
    t = solution.t
    p = plot(title="Phase evolution", xaxis="t [ms]")
    plot!(t, ψ_array[:, 1], label="ψ₁")
    plot!(t, ψ_array[:, 2], label="ψ₂")
    display(p)
    return nothing
end
