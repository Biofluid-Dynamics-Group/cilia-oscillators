using LinearAlgebra
using LinearSolve
using DifferentialEquations
using ODEInterfaceDiffEq
using Plots
using Base.Threads


include("filament_oscillator.jl")
include("../utils/greens_functions.jl")


M = 2
# φ = 12.5*π/180.0
φ = 0.0
A = [cos(φ) -sin(φ) 0.0; sin(φ) cos(φ) 0.0; 0.0 0.0 1.0]
d = 40.0

N = 5
h = L/(N - 1)
s = [(i)*h for i=1:N]
x₀ = [[(j - 1)*d, 0.0, 0.0] for j=1:M]

"""
    x(Ψ::Vector)

Returns the stacked positions of the discretised filaments at phase `Ψ`. It is an NM x 3
size matrix with the positions of the `j`th filament at the `i`th discretisation point
given by `x(Ψ)[i + (j - 1)N, :]`.
"""
function x(Ψ::Vector)
    result = zeros(N*M, 3)
    @threads for idx=1:M*N
        j = ceil(Int, idx / N)
        i = idx - (j - 1)*N
        result[idx, :] .= A*ξ(s[i], Ψ[2j - 1:2j]) + x₀[j]
    end
    return result
end

"""
    Ψ_dot(Ψ::Vector, μ::Real, a::Real, κ_b::Real)

Returns the phase velocity at a given system phase `Ψ`, in a fluid with viscosity `mu`,
where all filaments were discretised using spheres of radius `a`, and they share bending
stiffness `κ_b`. The phase velocity is determined by solving the saddle-point system that
serves as the equation of motion for the system of filaments.
"""
function Ψ_dot(Ψ::Vector, μ::Real, a::Real, κ_b::Real)
    Q = Q_ref(Ψ, μ, a, κ_b)
    κ_matrix = κₕ(Ψ)
    matrix = [Πₕ(Ψ, μ, a) -κ_matrix; -κ_matrix' zeros(2M, 2M)]
    rhs = vcat(zeros(M*3N), -Q)
    problem = LinearProblem(matrix, rhs)
    solution = solve(problem, KrylovJL_GMRES())
    dΨ_dt = solution[end-(2*M - 1):end]
    return dΨ_dt
end

"""
    block_diagonal(matrices::Array{Matrix})

Given an array of equally-sized matrices, returns a block-diagonal matrix with the
matrices of the array on the diagonal.
"""
function block_diagonal(matrices::Vector)
    example = matrices[1]
    rows, columns = size(example)
    n = length(matrices)
    result = zeros(rows*n, columns*n)
    for i = 1:n
        result[1 + (i - 1)*rows:i*rows, 1 + (i - 1)*columns:i*columns] .= matrices[i]
    end
    return result
end

"""
    κₕ(Ψ::Vector)

Returns the block-diagonal matrix that maps the phase velocity to the discretised filament
system velocity at phase `Ψ`.
"""
function κₕ(Ψ::Vector)
    return block_diagonal([Kₕ(Ψ[2j - 1:2j]) for j=1:M])
end

"""
    Πₕ(Ψ::Vector, μ::Real, a::Real)

Returns the block-diagonal regularised Stokes mobility matrix mapping forces in the
discretised filament system to velocities at phase `Ψ`, in a fluid with viscosity `μ`,
where the filaments were discretised using spheres of radius `a`.
"""
function Πₕ(Ψ::Vector, μ::Real, a::Real)
    mobility = zeros(3*N*M, 3*N*M)
    x_vector = x(Ψ)
    @threads for i = 1:N*M
        @threads for j = 1:N*M
            α = 1 + 3*(i - 1)
            β = 1 + 3*(j - 1)
            tensor = rotne_prager_blake_tensor(x_vector[i, :], x_vector[j, :], μ, a)
            if any(isnan, x_vector)
                @warn "NaN detected in x_vector"
            end    
            if any(isnan, tensor)
                @warn "NaN detected in rotne_prager_blake_tensor at i=$i, j=$j"
            end
            mobility[α:(α + 2), β:(β + 2)] .= tensor
        end
    end
    return mobility
end

"""
    Kₕ(ψ::Vector)

Returns the matrix that maps the phase velocity to the discretised filament velocity at
phase `ψ`.
"""
function Kₕ(ψ::Vector)
    return cat([A*K(s[i], ψ) for i = 1:N]..., dims=1)
end

"""
    Q_ref(Ψ::Vector, μ::Real, a::Real, κ_b::Real)

Returns the generalised force vector of forcings that, in absence of other cilia, induce a
reference beat in each filament. It includes the elastic relaxation.
"""
function Q_ref(Ψ::Vector, μ::Real, a::Real, κ_b::Real)
    function q_ref_indexed(j::Int)
        K = Kₕ(Ψ[2j - 1:2j])
        M = Mₕ(Ψ[2j - 1:2j], μ, a)
        return K'*(M\(K*[2π/T, 0.0])) .- [0.0, -κ_b*Ψ[2j]]
    end
    return cat([q_ref_indexed(j) for j=1:M]..., dims=1)
end

"""
    filament_oscillators!(dΨ_dt::Vector, Ψ::Vector, p::NamedTuple, t::Real)

Updates the phase velocities `dΨ_dt` given the current phase `Ψ` and the parameters `p`.
"""
function filament_oscillators!(dΨ_dt::Vector, Ψ::Vector, p::NamedTuple, t::Real)
    μ = p.μ
    a = p.a
    κ_b = p.κ_b
    rhs = Ψ_dot(Ψ, μ, a, κ_b)
    dΨ_dt .= rhs
    return nothing
end

"""
    run_system(p::NamedTuple, final_time::Real, num_steps::Int)

Runs the system of filaments with parameters `p` for `final_time` seconds, using `num_steps` timesteps with 4th order Adams-Bashforth.
"""
function run_system(p::NamedTuple, final_time::Real, num_steps::Int, alg)
    t_span = (0.0, final_time)
    Ψ₀ = repeat([1.0, 0.0]; outer=[M])
    problem = ODEProblem(filament_oscillators!, Ψ₀, t_span, p)
    solution = solve(problem, alg, dt=final_time/num_steps)
    return solution
end

