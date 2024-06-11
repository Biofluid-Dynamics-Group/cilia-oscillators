include("../src/utils/integral.jl")
using .Integral
using ProgressMeter
using Plots

include("../src/beats/tangent_angle.jl")


L = 20.0                  # Cilium length
T = 60.0                  # Beat period
θ_0 = π/2.1               # Maximum tangent angle
f_eff = 0.3               # Effective stroke fraction
f_ψ = 0.85                # Fraction of recovery controlled by travelling wave
f_w = 0.4                 # Travelling wavelength w/r to cilium length
orientation = π/2.0       # Cilium orientation
rule = GaussLegendre(15)  # Quadrature rule for numerical integration
beat_params = BeatParameters(L, T, θ_0, f_eff, f_ψ, f_w, orientation, rule)
ψ₁_array = collect(0.0:0.1:2π)
ψ₂_array = collect(-π:0.1:π)
s_array = collect(0.0:0.1:L)

color_scheme_1 = palette(:twilight, length(ψ₁_array), rev=true)
color_scheme_2 = palette(:twilight, length(ψ₂_array), rev=true)

p1 = plot(xlims=(-L, L), ylims=(0.0, L), legend=false)
@gif for (k, ψ₁) in enumerate(ψ₁_array)
    positions = zeros(length(s_array), 3)
    for (i, s) in enumerate(s_array)
        positions[i, :] .= ξ(s, [ψ₁, 0.0], beat_params)
    end
    plot!(
        p1, positions[:, 1], positions[:, 3], title="ψ₁ variation", color=color_scheme_1[k]
    )
end

p2 = plot(xlims=(-L, 20.0), ylims=(-L, L), legend=false)
@gif for (k, ψ₂) in enumerate(ψ₂_array)
    positions = zeros(length(s_array), 3)
    for (i, s) in enumerate(s_array)
        positions[i, :] .= ξ(s, [π*f_eff, ψ₂], beat_params)
    end
    plot!(p2, positions[:, 1], positions[:, 3], color=color_scheme_2[k])
end

p = plot(p1, p2, layout=(1, 2), title="ψ₂ variation")
display(p)
