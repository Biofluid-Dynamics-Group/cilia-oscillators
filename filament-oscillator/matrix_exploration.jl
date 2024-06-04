include("src/utils/integral.jl")
using .Integral
using ProgressMeter

include("src/beats/tangent_angle.jl")
include("src/models/physical_params.jl")


# Integration rule
rule = GaussLegendre(15)

# Instance objects
beat_params = BeatParameters(L, T, θ_0, f_eff, f_ψ, f_w, orientation, rule)

# ψ₁_array = collect(0.0:0.1:2π)
# ψ₂_array = collect(-π:0.1:π)
s = collect(0.0:0.1:L)
# matrix_array = zeros(length(ψ₁_array), length(ψ₂_array), 2, 2)
# @showprogress for (i, ψ₁) in enumerate(ψ₁_array)
#     for (j, ψ₂) in enumerate(ψ₂_array)
#         K_matrix = vcat([K(s[i], [ψ₁, ψ₂], beat_params) for i=1:length(s)]...)
#         matrix_array[i, j, :, :] = K_matrix'*K_matrix
#     end
# end

function det_of_KTK(ψ₁, ψ₂)
    K_matrix = vcat([K(s[i], [ψ₁, ψ₂], beat_params) for i=1:length(s)]...)
    return det(K_matrix'*K_matrix)/25000
end

using Plots

xs = range(0.0, stop=2π, length=100)
ys = range(-π, stop=π, length=100)
contourf(xs, ys, det_of_KTK, xlabel="ψ₁", ylabel="ψ₂", zlabel="det(KᵀK)")

# plot(xs, det_of_KTK.(xs, 0.0), xlabel="ψ₁", ylabel="det(KᵀK)", label="ψ₂=0.0")
