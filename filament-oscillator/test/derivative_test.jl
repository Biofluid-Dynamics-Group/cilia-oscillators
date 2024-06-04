include("../src/utils/integral.jl")
using .Integral
using ProgressMeter
using Test

include("../src/beats/tangent_angle.jl")
include("../src/models/physical_params.jl")

rule = GaussLegendre(140)
params = BeatParameters(L, T, θ_0, f_eff, f_ψ, f_w, orientation, rule)

@testset "Derivative of θ_rec" begin
    for s = 0.0:1:L
        ∂θ_at_s(ψ) = ∂θ_rec_∂ψ_1(s, ψ, params)
        for ψ₁ = 0.0:0.5:2π*params.f_rec
            @test ∫(0.0, ψ₁, ∂θ_at_s, rule) ≈ θ_rec(s, ψ₁, params) - θ_rec(s, 0.0, params) atol=1e-6
        end
    end
end

@testset "Derivative of θ_eff" begin
    for ψ₁ = 0.0:0.1:2π*params.f_eff
        ∂θ_no_params(ψ) = ∂θ_eff_∂ψ_1(ψ, params)
        @test ∫(0.0, ψ₁, ∂θ_no_params, rule) ≈ θ_eff(ψ₁, params) - θ_eff(0.0, params) atol=1e-6
    end
end
