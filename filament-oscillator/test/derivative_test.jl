using Test

include("../src/beats/tangent_angle.jl")

rule = GaussLegendre(140)

@testset "Derivative of θ_rec" begin
    for s = 0.0:1:L
        ∂θ_at_s(ψ) = ∂θ_rec_∂ψ_1(s, ψ)
        for ψ₁ = 0.0:0.5:2π*f_rec
            @test ∫(0.0, ψ₁, ∂θ_at_s, rule) ≈ θ_rec(s, ψ₁) - θ_rec(s, 0.0) atol=1e-6
        end
    end
end

@testset "Derivative of θ_eff" begin
    for ψ₁ = 0.0:0.1:2π*f_eff
        @test ∫(0.0, ψ₁, ∂θ_eff_∂ψ_1, rule) ≈ θ_eff(ψ₁) - θ_eff(0.0) atol=1e-6
    end
end
