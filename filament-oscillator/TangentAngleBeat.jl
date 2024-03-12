using LinearAlgebra
using FastGaussQuadrature


order = 20  # Quadrature order
nodes, weights = gausslegendre(order)

T = 60.
L = 19.
θ_0 = π/2.1
f_eff = 0.3
f_ψ = 0.85
f_w = 0.4
orientation = π/2.0

T_eff = T*f_eff
T_rec = T - T_eff
w = f_w*L
c = (L + w)/T_rec

function θ(s::Real, t::Real)
    period_t = mod(t, T)
    if period_t < T_eff
        return θ_eff(t)
    else
        return θ_rec(s, period_t - T_eff)
    end
end

function θ_eff(t::Real)
    return θ_0*cos(π*t/T_eff)
end

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

function θ_rec(s::Real, t::Real)
    return θ_0*((1.0 - f_ψ)*sin(π*t/T_rec - π/2.0) - f_ψ*g(
        (s - c*t)/w + 0.5
    ))
end

function ∫(a::Real, b::Real, f::Function, method::String)
    if b - a < eps()
        return 0.0
    end
    if method == "gauss"
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

function tangent_at_t(t::Real)
    function tan_x(s::Real)
        return cos(θ(s, t) + orientation)
    end
    function tan_z(s::Real)
        return sin(θ(s, t) + orientation)    
    end
    return tan_x, tan_z
end

function ξ_time(s::Real, t::Real)
    tan_x, tan_z = tangent_at_t(t)
    if s < eps()
        return [0., 0., 0.]
    end
    x_position = ∫(0., s, tan_x, "trapezoidal")
    z_position = ∫(0., s, tan_z, "trapezoidal")
    return [x_position, 0., z_position]
end

function ξ_time(s::Vector, t::Real)
    tan_x, tan_z = tangent_at_t(t)
    N, = size(s)
    return L/N*[cumsum(tan_x.(s)), 0., cumsum(tan_z.(s))]
end

function ξ(s::Real, ψ_1::Real)
    t = T*ψ_1/(2.0*π)
    return ξ_time(s, t)
end

function ξ(s::Vector, ψ_1::Real)
    t = T*ψ_1/(2.0*π)
    return ξ_time(s, t)
end

function dθ_dt(s::Real, t::Real)
    period_t = mod(t, T)
    if period_t < T_eff
        return dθ_eff_dt(t)
    else
        return dθ_rec_dt(s, t)
    end
end

function dθ_eff_dt(t::Real)
    return -θ_0*sin(π*t/T_eff)*π/T_eff
end

function dg_dx(u::Real)
    if u ≤ 0.5 || u ≥ 0.5
        return 0.0
    else
        return -4.0*(u^2 + 1.0)*sech(4.0*u/(4.0*u^2 - 1.0))^2/((1.0 - 4.0*u^2)^2)
    end
end

function dθ_rec_dt(s::Real, t::Real)
    first_term = (1.0 - f_ψ)*cos(π*t/T_rec + π/2.0)*π/T_rec
    second_term = f_ψ*dg_dx((s - c*t)/w + 1.0/2.0)*c/w
    return θ_0*(first_term + second_term)
end

function normal_at_t(t::Real)
    function n_x(s::Real)
        return -sin(θ(s, t) + orientation)
    end
    function n_z(s::Real)
        return cos(θ(s, t) + orientation)
    end
    return n_x, n_z
end

function dθ_dt_at_t(t::Real)
    function dθ_dt_arclength(s::Real)
        return dθ_dt(s, t)
    end
    return dθ_dt_arclength
end

function dξ_dψ_1_time(s::Real, t::Real)
    jacobian = T/(2.0*π)
    dθ = dθ_dt_at_t(t)
    n_x, n_z = normal_at_t(t)
    function integrand_x(s_prime::Real)
        return n_x(s_prime)*dθ(s_prime)*jacobian
    end
    function integrand_z(s_prime::Real)
        return n_z(s_prime)*dθ(s_prime)*jacobian
    end
    x_position = ∫(0., s, integrand_x, "trapezoidal")
    z_position = ∫(0., s, integrand_z, "trapezoidal")
    return [x_position, 0., z_position]
end

function dξ_dψ_1(s::Real, ψ_1::Real)
    t = T*ψ_1/(2.0*π)
    return dξ_dψ_1_time(s, t)
end
