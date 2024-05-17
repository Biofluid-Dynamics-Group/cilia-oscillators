using Plots

include("../src/beats/tangent_angle.jl")


ψ₁_array = collect(0.0:0.1:4π)
ψ₂_array = collect(-π:0.1:π)
s_array = collect(0.0:0.1:20.0)

p = plot(xlims=(-20.0, 20.0), ylims=(0.0, 20.0), legend=false)
@gif for ψ₁ in ψ₁_array
    positions = zeros(length(s_array), 3)
    for (i, s) in enumerate(s_array)
        positions[i, :] .= ξ(s, [ψ₁, 0.0])
    end
    plot(p, positions[:, 1], positions[:, 3])
end
display(p)

p = plot(xlims=(-20.0, 20.0), ylims=(-20.0, 20.0), legend=false)
@gif for ψ₂ in ψ₂_array
    positions = zeros(length(s_array), 3)
    for (i, s) in enumerate(s_array)
        positions[i, :] .= ξ(s, [π*f_eff, ψ₂])
    end
    plot(p, positions[:, 1], positions[:, 3])
end
display(p)
