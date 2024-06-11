include("../src/utils/integral.jl")
using .Integral
using ProgressMeter
using Test

include("../src/beats/tangent_angle.jl")


L = 20.0                  # Cilium length
T = 60.0                  # Beat period
θ_0 = π/2.1               # Maximum tangent angle
f_eff = 0.3               # Effective stroke fraction
f_ψ = 0.85                # Fraction of recovery controlled by travelling wave
f_w = 0.4                 # Travelling wavelength w/r to cilium length
orientation = π/2.0       # Cilium orientation
rule = GaussLegendre(20)  # Quadrature rule for numerical integration
params = BeatParameters(L, T, θ_0, f_eff, f_ψ, f_w, orientation, rule)

atol = 1e-8

@testset "Derivative of θ_rec" begin
    p = Progress(ceil(Int, L÷1 + 2π*params.f_rec÷0.5); desc="Testing ∂θ/∂ψ₁ in recovery stroke...")
    for s = 0.0:1:L
        ∂θ_at_s(ψ) = ∂θ_rec_∂ψ_1(s, ψ, params)
        for ψ₁ = 0.0:0.5:2π*params.f_rec
            @test ∫(0.0, ψ₁, ∂θ_at_s, rule) ≈ θ_rec(s, ψ₁, params) - θ_rec(s, 0.0, params) atol=atol
            next!(p)
        end
    end
end

@testset "Derivative of θ_eff" begin
    p = Progress(ceil(Int, 2π*params.f_eff÷0.1); desc="Testing ∂θ/∂ψ₁ in effective stroke...")
    for ψ₁ = 0.0:0.1:2π*params.f_eff
        ∂θ_no_params(ψ) = ∂θ_eff_∂ψ_1(ψ, params)
        @test ∫(0.0, ψ₁, ∂θ_no_params, rule) ≈ θ_eff(ψ₁, params) - θ_eff(0.0, params) atol=atol
        next!(p)
    end
end

