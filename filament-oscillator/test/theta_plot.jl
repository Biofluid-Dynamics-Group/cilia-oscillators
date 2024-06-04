include("../src/utils/integral.jl")
using .Integral
using ProgressMeter
using Test

include("../src/beats/tangent_angle.jl")
include("../src/models/physical_params.jl")

rule = GaussLegendre(140)
params = BeatParameters(L, T, θ_0, f_eff, f_ψ, f_w, orientation, rule)

ψ₁ = 0.0:0.01:4π
# s = 0.0:0.1:L

p = plot()

# for s_val in s
∂θ_plot(ψ) = ∂θ_∂ψ_1(0.0, ψ, params)
plot!(ψ₁, ∂θ_plot.(ψ₁), xlabel="ψ₁", ylabel="∂θ/∂ψ₁")
# end

display(p)
