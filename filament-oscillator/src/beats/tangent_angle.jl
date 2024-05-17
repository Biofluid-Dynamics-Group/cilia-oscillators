using LinearAlgebra

include("../utils/integral.jl")
include("../models/physical_params.jl")


rule = Trapezoidal(20)

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
        return cos(θ(σ, ψ[1]) + ψ[2]*σ/L + orientation)
    end
    function z_integrand(σ::Real)
        return sin(θ(σ, ψ[1]) + ψ[2]*σ/L + orientation)
    end
    x = ∫(0.0, s, x_integrand, rule)
    z = ∫(0.0, s, z_integrand, rule)
    return [x, 0.0, z]
end

"""
    ∂ξ_∂ψ_1(s::Real, ψ::Vector)

Returns the derivative of the position of the filament at arglenth `s` and phase `ψ` with
respect to `ψ_1`.
"""
function ∂ξ_∂ψ_1(s::Real, ψ::Vector)
    function x_integrand(σ::Real)
        return -sin(θ(σ, ψ[1]) + ψ[2]*σ/L + orientation)*∂θ_∂ψ_1(σ, ψ[1])
    end
    function z_integrand(σ::Real)
        return cos(θ(σ, ψ[1]) + ψ[2]*σ/L + orientation)*∂θ_∂ψ_1(σ, ψ[1])
    end
    x = ∫(0.0, s, x_integrand, rule)
    z = ∫(0.0, s, z_integrand, rule)
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
    x = ∫(0.0, s, x_integrand, rule)
    z = ∫(0.0, s, z_integrand, rule)
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
