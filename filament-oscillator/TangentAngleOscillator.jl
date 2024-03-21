using LinearAlgebra
using FastGaussQuadrature
using Plots


order = 10  # Quadrature order
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

function θ_bend(s::Real, ψ_1::Real)
    k = ψ_1 ÷ (2.0*π)
    if mod(ψ_1, 2π) < f_eff*2.0*π
        return (2*k - 2)*θ_0*f_ψ
    else
        t = mod(ψ_1, 2π)*T/(2.0*π) - T_eff
        return (2*k - 1)*θ_0*f_ψ - θ_0*f_ψ*g((s - c*t)/w + 0.5)
    end
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

function tangent_at_ψ_1(ψ_1::Real)
    function tan_x(s::Real)
        return cos(θ_bend(s, ψ_1))
    end
    function tan_z(s::Real)
        return sin(θ_bend(s, ψ_1))    
    end
    return tan_x, tan_z
end

function normal_at_ψ_1(ψ_1::Real)
    function nor_x(s::Real)
        return -sin(θ_bend(s, ψ_1))
    end
    function nor_z(s::Real)
        return cos(θ_bend(s, ψ_1))
    end
    return nor_x, nor_z
end

function ξ(s::Real, ψ_1::Real)
    tan_x, tan_z = tangent_at_ψ_1(ψ_1)
    x_position = ∫(0.0, s, tan_x, "trapezoidal")
    z_position = ∫(0.0, s, tan_z, "trapezoidal")
    return [x_position, 0.0, z_position]
end

function dθ_bend_dψ_1(s::Real, ψ_1::Real)
    k = ψ_1 ÷ (2.0*π)
    if mod(ψ_1, 2.0*π) < f_eff*2.0*π
        return 0.0
    else
        t = mod(ψ_1, 2π)*T/(2.0*π) - T_eff
        return -(k - 1)*θ_0*f_ψ*dg_dx((s - c*t)/w + 0.5)*c*T/(2.0*π*w)
    end
end

function dg_dx(u::Real)
    if u ≤ -0.5 || u ≥ 0.5
        return 0.0
    else
        return 4.0*(4.0*u^2 + 1.0)*sech(4.0*u/(4.0*u^2 - 1.0))^2/((1.0 - 4.0*u^2)^2)
    end
end

function dξ_dψ_1(s::Real, ψ_1::Real)
    nor_x, nor_z = normal_at_ψ_1(ψ_1)
    function x_integrand(s::Real)
        return nor_x(s)*dθ_bend_dψ_1(s, ψ_1)
    end
    function z_integrand(s::Real)
        return nor_z(s)*dθ_bend_dψ_1(s, ψ_1)
    end
    x_position = ∫(0.0, s, x_integrand, "trapezoidal")
    z_position = ∫(0.0, s, z_integrand, "trapezoidal")
    return [x_position, 0.0, z_position]
end

# p = plot(title="θ_bend")
# # for s ∈ range(0.0, L, 1)
# s = L/2
# θ_arr = θ_bend.(s, 0:0.1:10*2π)
# ψ_arr = collect(0:0.1:10*2π)
# plot!(p, ψ_arr, θ_arr, label="s = $s")
# # end
# display(p)

p = plot(title="θ_bend in s")
for ψ_1 ∈ range(0.0, 2*pi, 10)
    θ_arr = θ_bend.(collect(range(0.0, L, 100)), ψ_1)
    plot!(p, collect(range(0.0, L, 100)), θ_arr, label="ψ_1 = $ψ_1")
end

display(p)