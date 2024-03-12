using LinearAlgebra
using LinearSolve
using DifferentialEquations
using Plots


include("TangentAngleBeat.jl")
include("../stokes/GreensFunctions.jl")

μ = 1.0
a = 7e-2

N = 5
h = L/(N - 1)
s = [(i)*h for i=1:N]

function R(ψ_2::Real)
    return [cos(ψ_2) 0 -sin(ψ_2); 0 1 0; sin(ψ_2) 0 cos(ψ_2)]
end

function x(s::Real, ψ::Vector)
    return R(ψ[2])*ξ(s, ψ[1])
end

function dR_dψ_2(ψ_2::Real)
    return [-sin(ψ_2) 0 -cos(ψ_2); 0 1 0; cos(ψ_2) 0 -sin(ψ_2)]
end

function K(s::Real, ψ::Vector)
    return [R(ψ[2])*dξ_dψ_1(s, ψ[1]) dR_dψ_2(ψ[2])*ξ(s, ψ[1])]
end

function Kₕ(ψ::Vector)
    return cat([K(s[i], ψ) for i = 1:N]..., dims=1)
end

function Mₕ(ψ::Vector, μ::Real, a::Real)
    mobility = zeros(3*N, 3*N)
    for i = 1:N
        for j = 1:N
            α = 1 + 3*(i - 1)
            β = 1 + 3*(j - 1)
            mobility[α:(α + 2), β:(β + 2)] .= rotne_prager_blake_tensor(
                x(s[i], ψ), x(s[j], ψ), μ, a
            )
        end
    end
    return mobility
end

function q_ref(ψ::Vector, μ::Real, a::Real)
    return Kₕ(ψ)'*inv(Mₕ(ψ, μ, a))*Kₕ(ψ)*[2*π/T, 0.0]
end

function ψ_dot(ψ::Vector, μ::Real, a::Real)
    M = Mₕ(ψ, μ, a)
    K = Kₕ(ψ)
    q = q_ref(ψ, μ, a)
    matrix = [M -K; -K' zeros(2, 2)]
    rhs = vcat(zeros(3*N), -q)
    problem = LinearProblem(matrix, rhs)
    solution = solve(problem, KrylovJL_GMRES())
    dψ_dt = solution.u[end-1:end]
    return dψ_dt
end

function filament_oscillator!(dψ_dt::Vector, ψ::Vector, p::NamedTuple, t::Real)
    μ = p.μ
    a = p.a
    rhs = ψ_dot(ψ, μ, a)
    dψ_dt[1] = rhs[1]
    dψ_dt[2] = rhs[2]
    return nothing
end

function run_filament(p::NamedTuple, final_time::Real, num_steps::Int)
    t_span = (0.0, final_time)
    problem = ODEProblem(filament_oscillator!, [2*π, 0.0], t_span, p)
    solution = solve(problem, AB4(), dt=final_time/num_steps)
    return solution
end

function plot_trajectory(solution, p::NamedTuple)
    ϕ = solution.u
    positions = zeros(size(ϕ)[1], 3)
    for i = 1:size(ϕ)[1]
        positions[i, :] .= x(ϕ[i], p.x_0)
    end
    display(plot!(positions[:, 1], positions[:, 3]))
    return nothing
end


p = (a=a, μ=μ)

solution = run_filament(p, 60.0, 100)
# plot_trajectory(solution, p)
