include("../src/utils/integral.jl")
using .Integral
using ProgressMeter
using Test

include("../src/beats/tangent_angle.jl")


L = 20.0                   # Cilium length
T = 60.0                   # Beat period
θ_0 = π/2.1                # Maximum tangent angle
f_eff = 0.3                # Effective stroke fraction
f_ψ = 0.85                 # Fraction of recovery controlled by travelling wave
f_w = 0.4                  # Travelling wavelength w/r to cilium length
orientation = π/2.0        # Cilium orientation
rule = Trapezoidal(25000)  # Quadrature rule for numerical integration
params = BeatParameters(L, T, θ_0, f_eff, f_ψ, f_w, orientation, rule)

atol = 1e-8

@testset "Derivative of ξ with respect to ψ₁" begin
    p = Progress(ceil(Int, 2π÷1 + L÷5 + 2π÷0.5); desc="Testing ∂ξ/∂ψ₁...")
    for ψ₂ = -π:1:π
        for s = 0.0:5:L
            ∂ξ1_at_s(ψ) = ∂ξ_∂ψ_1(s, [ψ, ψ₂], params)[1]
            ∂ξ2_at_s(ψ) = ∂ξ_∂ψ_1(s, [ψ, ψ₂], params)[2]
            ∂ξ3_at_s(ψ) = ∂ξ_∂ψ_1(s, [ψ, ψ₂], params)[3]
            for ψ₁ = 0.0:0.5:2π
                @test ∫(0.0, ψ₁, ∂ξ1_at_s, rule) ≈ ξ(s, [ψ₁, ψ₂], params)[1] - ξ(s, [0.0, ψ₂], params)[1] atol=atol
                @test ∫(0.0, ψ₁, ∂ξ2_at_s, rule) ≈ ξ(s, [ψ₁, ψ₂], params)[2] - ξ(s, [0.0, ψ₂], params)[2] atol=atol
                @test ∫(0.0, ψ₁, ∂ξ3_at_s, rule) ≈ ξ(s, [ψ₁, ψ₂], params)[3] - ξ(s, [0.0, ψ₂], params)[3] atol=atol
                next!(p)
            end
        end
    end
end

@testset "Derivative of ξ with respect to ψ₁" begin
    p = Progress(ceil(Int, 2π÷1 + L÷5 + 2π÷0.5); desc="Testing ∂ξ/∂ψ₂...")
    for ψ₁ = 0.0:1:2π
        for s = 0.0:5:L
            ∂ξ1_at_s(ψ) = ∂ξ_∂ψ_2(s, [ψ₁, ψ], params)[1]
            ∂ξ2_at_s(ψ) = ∂ξ_∂ψ_2(s, [ψ₁, ψ], params)[2]
            ∂ξ3_at_s(ψ) = ∂ξ_∂ψ_2(s, [ψ₁, ψ], params)[3]
            for ψ₂ = -π:0.5:π
                @test ∫(0.0, ψ₂, ∂ξ1_at_s, rule) ≈ ξ(s, [ψ₁, ψ₂], params)[1] - ξ(s, [0.0, ψ₂], params)[1] atol=atol
                @test ∫(0.0, ψ₂, ∂ξ2_at_s, rule) ≈ ξ(s, [ψ₁, ψ₂], params)[2] - ξ(s, [0.0, ψ₂], params)[2] atol=atol
                @test ∫(0.0, ψ₂, ∂ξ3_at_s, rule) ≈ ξ(s, [ψ₁, ψ₂], params)[3] - ξ(s, [0.0, ψ₂], params)[3] atol=atol
                next!(p)
            end
        end
    end
end
