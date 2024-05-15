using LinearAlgebra
using FastGaussQuadrature

include("../models/physical_params.jl")


const order = 10  # Quadrature order
nodes, weights = gausslegendre(order)

"""
    g(u::Real)

Returns the value of the transition function `g` for a given argument `u`.
"""
function g(u::Real)
    if u ≤ -0.5
        return -1.0
    elseif -0.5 < u && u < 0.5
        return 2.0*exp(-2.0/(2.0*u + 1.0))/(
            exp(-2.0/(2.0*u + 1.0)) + exp(2.0/(2.0*u - 1.0))
        ) - 1.0
    else
        return 1.0
    end
end

"""
    dg_dx(u::Real)

Returns the derivative of the transition function `g` at argument `u`.
"""
function dg_dx(u::Real)
    if u ≤ -0.5 || u ≥ 0.5
        return 0.0
    else
        return 4.0*(4.0*u^2 + 1.0)*sech(4.0*u/(4.0*u^2 - 1.0))^2/((1.0 - 4.0*u^2)^2)
    end
end

"""
    ∫(a::Real, b::Real, f::Function, method::String)

Numerically integrates the function `f` on the domain (`a`, `b`) with method `method`.
"""
function ∫(a::Real, b::Real, f::Function, method::String)
    if b - a < eps()
        return 0.0
    end
    if method == "gausslegendre"
        function fˢ(v::Real)
            u = (b - a)*v/2.0 + (a + b)/2.0
            return 2.0/(b - a)*f(u)
        end
        return dot(weights, fˢ.(nodes))
    elseif method == "trapezoidal"
        Δx = (b - a)/order
        x = [a + i*Δx for i=0:order]
        return Δx/2.0*dot(f.(x), vcat([1], 2*ones(order-1), [1]))
    end
end

"""
    θ(s::Real, ψ_1::Real)

Returns the preffered tangent angle for the given arclength `s` and shape phase `ψ_1`.
"""
function θ(s::Real, ψ_1::Real)
    mod_ψ_1 = mod2pi(ψ_1)
    if mod_ψ_1 < f_eff*2π
        return θ_eff(mod_ψ_1)
    else
        return θ_rec(s, mod_ψ_1 - f_eff*2π)
    end
end

"""
    θ_eff(ψ_1::Real)

Returns the angle of the filament during the effective stroke at shape phase `ψ_1`.
"""
function θ_eff(ψ_1::Real)
    return θ_0*cos(ψ_1/(2.0*f_eff))
end

"""
    θ_rec(s::Real, ψ_1::Real)

Returns the tangent angle of the filament during the recovery stroke at arclength `s` and
shape phase `ψ_1`.
"""
function θ_rec(s::Real, ψ_1::Real)
    return θ_0*((1.0 - f_ψ)*sin(ψ_1/(2.0*f_rec) - π/2.0) - f_ψ*g(
        (s - (c*T/2π)*ψ_1)/w + 0.5
    ))
end

"""
    ∂θ_∂ψ_1(s::Real, ψ_1::Real)

Returns the derivative of the tangent angle of the filament with respect to `ψ_1` at
arclength `s` and shape phase `ψ_1`.
"""
function ∂θ_∂ψ_1(s::Real, ψ_1::Real)
    mod_ψ_1 = mod2pi(ψ_1)
    if mod_ψ_1 < f_eff*2π
        value = ∂θ_eff_∂ψ_1(mod_ψ_1)
    else
        value = ∂θ_rec_∂ψ_1(s, mod_ψ_1 - f_eff*2π)
    end
    return value
end

"""
    ∂θ_eff_∂ψ_1(ψ_1::Real)

Returns the derivative of the tangent angle of the filament during the effective stroke
at shape phase ψ_1.
"""
function ∂θ_eff_∂ψ_1(ψ_1::Real)
    value = -θ_0*sin(ψ_1/(2.0*f_eff))/(2.0*f_eff) 
    if value < eps()
        println("∂θ/∂ψ₁ = 0")
    end
   return -θ_0*sin(ψ_1/(2.0*f_eff))/(2.0*f_eff) 
end

"""
    ∂θ_rec_∂ψ_1(s::Real, ψ_1::Real)

Returns the derivative of the tangent angle of the filament during the recovery stroke
with respect to `ψ_1` at arclength `s` and shape phase `ψ_1`.
"""
function ∂θ_rec_∂ψ_1(s::Real, ψ_1::Real)
    first_term = (1.0 - f_ψ)*cos(ψ_1/(2.0*f_rec) + π/2.0)/(2.0*f_rec)
    second_term = f_ψ*dg_dx((s - (c*T/2π)*ψ_1)/w + 1.0/2.0)*(c*T/2π)/w
    return θ_0*(first_term + second_term)
end

"""
    ξ(s::Real, ψ::Vector)

Returns the position of the filament at arglenth `s` and phase `ψ`.
"""
function ξ(s::Real, ψ::Vector)
    function x_integrand(σ::Real)
        return cos(θ(σ, ψ[1]) + ψ[2]*s/L + orientation)
    end
    function z_integrand(σ::Real)
        return sin(θ(σ, ψ[1]) + ψ[2]*s/L + orientation)
    end
    x = ∫(0.0, s, x_integrand, "trapezoidal")
    z = ∫(0.0, s, z_integrand, "trapezoidal")
    return [x, 0.0, z]
end

"""
    ∂ξ_∂ψ_1(s::Real, ψ::Vector)

Returns the derivative of the position of the filament at arglenth `s` and phase `ψ` with
respect to `ψ_1`.
"""
function ∂ξ_∂ψ_1(s::Real, ψ::Vector)
    function x_integrand(σ::Real)
        return -sin(θ(σ, ψ[1]) + ψ[2]*s/L + orientation)*∂θ_∂ψ_1(σ, ψ[1])
    end
    function z_integrand(σ::Real)
        return cos(θ(σ, ψ[1]) + ψ[2]*s/L + orientation)*∂θ_∂ψ_1(σ, ψ[1])
    end
    x = ∫(0.0, s, x_integrand, "trapezoidal")
    z = ∫(0.0, s, z_integrand, "trapezoidal")
    return [x, 0.0, z]
end

"""
    ∂ξ_∂ψ_2(s::Real, ψ::Vector)

Returns the derivative of the position of the filament at arglenth `s` and phase `ψ` with
respect to `ψ_2`.
"""
function ∂ξ_∂ψ_2(s::Real, ψ::Vector)
    function x_integrand(σ::Real)
        return -sin(θ(σ, ψ[1]) + ψ[2]*σ/L + orientation)*σ/L
    end
    function z_integrand(σ::Real)
        return cos(θ(σ, ψ[1]) + ψ[2]*σ/L + orientation)*σ/L
    end 
    x = ∫(0.0, s, x_integrand, "trapezoidal")
    z = ∫(0.0, s, z_integrand, "trapezoidal")
    return [x, 0.0, z]
end

"""
    K(s::Real, ψ::Vector)

Returns the matrix that maps the phase velocity to the filament velocity at arclength `s`
and phase `ψ` in the reference beat plane.
"""
function K(s::Real, ψ::Vector)
    return [∂ξ_∂ψ_1(s, ψ) ∂ξ_∂ψ_2(s, ψ)]
end
