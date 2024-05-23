module Integral


using FastGaussQuadrature
using LinearAlgebra

export QuadratureRule, Trapezoidal, GaussLegendre, ∫

abstract type QuadratureRule end

struct Trapezoidal <: QuadratureRule
    order::Int
end

struct GaussLegendre <: QuadratureRule
    order::Int
    weights::Vector{Real}
    nodes::Vector{Real}
    function GaussLegendre(order::Int)
        nodes, weights = gausslegendre(order)
        return new(order, weights, nodes)
    end
end


"""
    ∫(a::Real, b::Real, f::Function, rule::QuadratureRule)

Numerically integrates the function `f` on the domain (`a`, `b`) with rule `rule`.
"""
function ∫(a::Real, b::Real, f::Function, rule::Trapezoidal)
    if b - a < eps()
        return 0.0
    end
    Δx = (b - a)/rule.order
    x = [a + i*Δx for i=0:rule.order]
    return Δx/2.0*dot(f.(x), vcat([1], 2*ones(rule.order-1), [1]))
end

function ∫(a::Real, b::Real, f::Function, rule::GaussLegendre)
    if b - a < eps()
        return 0.0
    end
    function fˢ(v::Real)
        u = (b - a)*v/2.0 + (a + b)/2.0
        return (b - a)/2.0*f(u)
    end
    return dot(rule.weights, fˢ.(rule.nodes))
end


end
