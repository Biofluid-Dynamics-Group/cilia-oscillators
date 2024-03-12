using LinearAlgebra
using BlockArrays
using DifferentialEquations

include("GreensFunctions.jl")
include("DefaultParameters.jl")


function R_2(θ::Real)
    return [[cos(θ) 0.0 -sin(θ)]; [0.0 1.0 0.0]; [sin(θ) 0.0 cos(θ)]]
end

function plane_projection(u::Vector, n::Vector)
    return u - (dot(u, n)/dot(n, n))*n
end

function x(φ::Real, r::Real, ζ::Real, x_0::Vector)
    return x_0 + [r*sin(φ), ζ, r*cos(φ)]
end

function e_r(x::Vector, x_0=x_0::Vector)
    vector = plane_projection(x - x_0, [0.0; 1.0; 0.0])
    return vector/norm(vector)
end

function e_φ(x::Vector, x_0=x_0::Vector)
    vector = cross(e_r(x, x_0), e_ζ)
    return vector/norm(vector)
end

const e_ζ = [0, -1, 0]

function r(x::Vector, x_0=x_0::Vector)
    plane_vector = [(x - x_0)[1], (x - x_0)[3]]
    return norm(plane_vector)
end

function ζ(x::Vector, x_0=x_0::Vector)
    return (x - x_0)[2]
end

function gamma(x::Vector, gamma_0=gamma_0::Real, a=a::Real)
    return gamma_0*(I + ((9*a)/(16*x[3]))*(I + [0, 0, 1]*[0, 0, 1]'))
end

function rotor!(dx_dt::Vector, x::Vector, p::NamedTuple, t::Real)
    λ = p.λ
    η = p.λ
    f_drive = p.f_drive
    r_0 = p.r_0
    a = p.a
    x_0 = p.x_0
    gamma_0 = p.gamma_0
    forces = -λ*(r(x, x_0) - r_0)*e_r(x,  x_0)
    forces += -η*ζ(x, x_0)*e_ζ
    forces += f_drive*e_φ(x, x_0)
    rhs = gamma(x, gamma_0, a)\forces
    dx_dt[1] = rhs[1]
    dx_dt[2] = rhs[2]
    dx_dt[3] = rhs[3]
    return nothing
end

function simulate_rotor(p::NamedTuple, final_time::Real, num_steps::Int)
    r_0 = p.r_0
    x_0 = p.x_0
    initial_x = x_0 + r_0*[0, 0, 1]
    t_span = (0.0, final_time)
    f = ODEFunction(rotor!)
    problem = ODEProblem(f, initial_x, t_span, p)
    solution = solve(problem, AB4(), dt=final_time/num_steps)
    return solution
end
