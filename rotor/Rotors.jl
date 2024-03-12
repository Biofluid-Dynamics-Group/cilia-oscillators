using LinearAlgebra
using DifferentialEquations
using Plots


function x(ϕ::Vector, x_0::Vector)
    r, φ, ζ = ϕ
    return x_0 + [r*sin(φ), ζ, r*cos(φ)]
end

function e_r(ϕ::Vector)
    return [sin(ϕ[2]), 0., cos(ϕ[2])]
end

function e_φ(ϕ::Vector)
    return [cos(ϕ[2]), 0., -sin(ϕ[2])]
end

e_ζ = [0., -1., 0.]

function polar_to_cartesian(ϕ::Vector)
    hcat(e_r(ϕ), e_φ(ϕ), e_ζ)
end

function gamma(x::Vector, gamma_0=gamma_0::Real, a=a::Real)
    return gamma_0*(I + ((9.0*a)/(16.0*x[3]))*(I + [0., 0., 1.]*[0., 0., 1.]'))
end

function rotor!(dϕ_dt::Vector, ϕ::Vector, p::NamedTuple, t::Real)
    r, φ, ζ = ϕ
    F = [-p.λ*(r - p.r_0), p.f, -p.η*ζ]
    Q = polar_to_cartesian(ϕ)
    friction = Q'*gamma(x(ϕ, p.x_0), p.gamma_0, p.a)*Q
    rhs = friction\F
    dϕ_dt[1] = rhs[1]
    dϕ_dt[2] = rhs[2]/r
    dϕ_dt[3] = rhs[3]
    return nothing
end

function run_rotor(p::NamedTuple, t_span::Tuple, num_steps::Integer)
    f = ODEFunction(rotor!)
    problem = ODEProblem(f, [p.r_0, pi/4, 0.], t_span, p)
    solution = solve(problem, AB4(), dt=t_span[2]/num_steps)
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

Λ = 0.1
num_beats = 5
R = 0.5
L = 2.
A = 0.01

μ = 1e3
d = 7.
T = 1/33

r_0 = R*d
a = A*d
gamma_0 = 6*π*μ*a
f = 2*π*r_0*gamma_0/T
λ = Λ*f/d
η = λ
x_0 = [0., 0., d]

p = (r_0=r_0, a=a, gamma_0=gamma_0, f=f, λ=λ, η=η, x_0=x_0)

solution = run_rotor(p, (0., num_beats*T), 1000)
plot_trajectory(solution, p)
