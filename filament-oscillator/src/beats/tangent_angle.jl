using LinearAlgebra


struct BeatParameters
    L::Real               # Cilium length
    T::Real               # Beat period
    θ_0::Real             # Maximum tangent angle
    f_eff::Real           # Effective stroke fraction
    f_ψ::Real             # Fraction of recovery controlled by travelling wave
    f_w::Real             # Travelling wavelength w/r to cilium length
    orientation::Real     # Normal angle to the beat plane
    rule::QuadratureRule  # Quadrature rule for numerical integration
    T_eff::Real           # Effective stroke duration
    f_rec::Real           # Recovery stroke fraction
    T_rec::Real           # Recovery stroke duration
    w::Real               # Travelling wavelength
    c::Real               # Travelling wave speed
    ψ_eff::Real           # Effective stroke phase
    ψ_rec::Real           # Recovery stroke phase
    function BeatParameters(
        L::Real, T::Real, θ_0::Real, f_eff::Real, f_ψ::Real, f_w::Real,
        orientation::Real, rule::QuadratureRule
    )
        T_eff = T*f_eff
        f_rec = 1.0 - f_eff
        T_rec = T*f_rec
        w = f_w*L
        c = (L + w)/T_rec
        ψ_eff = 2π*f_eff
        ψ_rec = 2π*f_rec
        new(
            L, T, θ_0, f_eff, f_ψ, f_w, orientation, rule, T_eff, f_rec, T_rec, w, c,
            ψ_eff, ψ_rec
        )
    end
end

"""
    g(x::Real)

Returns the value of the transition function `g` for a given argument `x`.
"""
function g(x::Real)
    if x ≤ -0.5
        return -1.0
    elseif -0.5 < x && x < 0.5
        return 2.0*exp(-2.0/(2.0*x + 1.0))/(
            exp(-2.0/(2.0*x + 1.0)) + exp(2.0/(2.0*x - 1.0))
        ) - 1.0
    else
        return 1.0
    end
end

"""
    dg_dx(x::Real)

Returns the derivative of the transition function `g` at argument `x`.
"""
function dg_dx(x::Real)
    if x ≤ -0.5 || x ≥ 0.5
        return 0.0
    else
        return 4.0*(4.0*x^2 + 1.0)*sech(4.0*x/(4.0*x^2 - 1.0))^2/((1.0 - 4.0*x^2)^2)
    end
end

"""
    θ(s::Real, ψ_1::Real)

Returns the preffered tangent angle for the given arclength `s` and shape phase `ψ_1`.
"""
function θ(s::Real, ψ_1::Real, params::BeatParameters)
    mod_ψ_1 = mod2pi(ψ_1)
    if mod_ψ_1 < params.ψ_eff
        return θ_eff(mod_ψ_1, params)
    else
        return θ_rec(s, mod_ψ_1 - params.ψ_eff, params)
    end
end

"""
    θ_eff(ψ_1::Real)

Returns the angle of the filament during the effective stroke at shape phase `ψ_1`.
"""
function θ_eff(ψ_1::Real, params::BeatParameters)
    # return params.θ_0*cos(0.5*ψ_1/params.f_eff)
    2params.θ_0/(2π*params.f_eff)*ψ_1
end

"""
    θ_rec(s::Real, ψ_1::Real)

Returns the tangent angle of the filament during the recovery stroke at arclength `s` and
shape phase `ψ_1`.
"""
function θ_rec(s::Real, ψ_1::Real, params::BeatParameters)
    # return params.θ_0*((1.0 - params.f_ψ)*sin(0.5*ψ_1/params.f_rec - 0.5*π) - params.f_ψ*g(
    #     (s - (params.c*params.T*ψ_1/2π))/params.w + 0.5
    # ))
    # return -params.θ_0*cos(0.5*ψ_1/params.f_rec)
    -2params.θ_0/(2π*params.f_rec)*ψ_1
end

"""
    ∂θ_∂ψ_1(s::Real, ψ_1::Real)

Returns the derivative of the tangent angle of the filament with respect to `ψ_1` at
arclength `s` and shape phase `ψ_1`.
"""
function ∂θ_∂ψ_1(s::Real, ψ_1::Real, params::BeatParameters)
    mod_ψ_1 = mod2pi(ψ_1)
    if mod_ψ_1 < params.ψ_eff
        value = ∂θ_eff_∂ψ_1(mod_ψ_1, params)
    else
        value = ∂θ_rec_∂ψ_1(s, mod_ψ_1 - params.ψ_eff, params)
    end
    return value
end

"""
    ∂θ_eff_∂ψ_1(ψ_1::Real)

Returns the derivative of the tangent angle of the filament during the effective stroke
at shape phase ψ_1.
"""
function ∂θ_eff_∂ψ_1(ψ_1::Real, params::BeatParameters)
#    return -0.5*params.θ_0*sin(0.5*ψ_1/params.f_eff)/params.f_eff 
    return 2params.θ_0/(2π*params.f_eff)
end

"""
    ∂θ_rec_∂ψ_1(s::Real, ψ_1::Real)

Returns the derivative of the tangent angle of the filament during the recovery stroke
with respect to `ψ_1` at arclength `s` and shape phase `ψ_1`.
"""
function ∂θ_rec_∂ψ_1(s::Real, ψ_1::Real, params::BeatParameters)
    # first_term = (1.0 - f_ψ)*0.5*cos(0.5*ψ_1/params.f_rec - 0.5*π)/params.f_rec
    # second_term = params.f_ψ*dg_dx((s - (params.c*params.T*ψ_1/2π))/params.w + 0.5)*(
    #     -params.c*params.T/2π/params.w
    # )
    # return params.θ_0*(first_term - second_term)
    # return 0.5*params.θ_0*sin(0.5*ψ_1/params.f_rec)/params.f_rec 
    return -2params.θ_0/(2π*params.f_rec)
end

"""
    ξ(s::Real, ψ::Vector)

Returns the position of the filament at arglenth `s` and phase `ψ`.
"""
function ξ(s::Real, ψ::Vector, params::BeatParameters)
    function x_integrand(σ::Real)
        return cos(θ(σ, ψ[1], params) + ψ[2]*σ/params.L + params.orientation)
    end
    function z_integrand(σ::Real)
        return sin(θ(σ, ψ[1], params) + ψ[2]*σ/params.L + params.orientation)
    end
    x = ∫(0.0, s, x_integrand, params.rule)
    z = ∫(0.0, s, z_integrand, params.rule)
    return [x, 0.0, z]
end

"""
    ∂ξ_∂ψ_1(s::Real, ψ::Vector)

Returns the derivative of the position of the filament at arglenth `s` and phase `ψ` with
respect to `ψ_1`.
"""
function ∂ξ_∂ψ_1(s::Real, ψ::Vector, params::BeatParameters)
    function x_integrand(σ::Real)
        return -sin(
            θ(σ, ψ[1], params) + ψ[2]*σ/params.L + params.orientation
        )*∂θ_∂ψ_1(σ, ψ[1], params)
    end
    function z_integrand(σ::Real)
        return cos(
            θ(σ, ψ[1], params) + ψ[2]*σ/params.L + params.orientation
        )*∂θ_∂ψ_1(σ, ψ[1], params)
    end
    x = ∫(0.0, s, x_integrand, params.rule)
    z = ∫(0.0, s, z_integrand, params.rule)
    return [x, 0.0, z]
end

"""
    ∂ξ_∂ψ_2(s::Real, ψ::Vector)

Returns the derivative of the position of the filament at arglenth `s` and phase `ψ` with
respect to `ψ_2`.
"""
function ∂ξ_∂ψ_2(s::Real, ψ::Vector, params::BeatParameters)
    function x_integrand(σ::Real)
        return -sin(θ(σ, ψ[1], params) + ψ[2]*σ/params.L + params.orientation)*σ/params.L
    end
    function z_integrand(σ::Real)
        return cos(θ(σ, ψ[1], params) + ψ[2]*σ/params.L + params.orientation)*σ/params.L
    end 
    x = ∫(0.0, s, x_integrand, params.rule)
    z = ∫(0.0, s, z_integrand, params.rule)
    return [x, 0.0, z]
end

"""
    K(s::Real, ψ::Vector)

Returns the matrix that maps the phase velocity to the filament velocity at arclength `s`
and phase `ψ` in the reference beat plane.
"""
function K(s::Real, ψ::Vector, params::BeatParameters)
    return [∂ξ_∂ψ_1(s, ψ, params) ∂ξ_∂ψ_2(s, ψ, params)]
end
