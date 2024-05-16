using Test

include("../src/utils/integral.jl")

atol = 1e-3

@testset "Trapezoidal" begin
    trapezoidal = Trapezoidal(42)

    # Test polyonomial integration
    @test ∫(0.0, 1.0, x -> x^2, trapezoidal) ≈ 1/3 atol=atol
    
    # Test integration of constant function
    @test ∫(0.0, 1.0, x -> 2.0, trapezoidal) ≈ 2.0 atol=atol
    
    # Test integration of sine function
    @test ∫(0.0, π, sin, trapezoidal) ≈ 2.0 atol=atol
    
    # Test integration of exponential function
    @test ∫(0.0, 1.0, exp, trapezoidal) ≈ exp(1.0) - 1.0 atol=atol
end

@testset "Gauss-Legendre" begin
    gausslegendre = GaussLegendre(10)

    # Test default interval
    @test ∫(-1.0, 1.0, exp, gausslegendre) ≈ exp(1.0) - exp(-1.0) atol=atol

    # Test polynomial integration
    @test ∫(0.0, 1.0, x -> x^2, gausslegendre) ≈ 1/3 atol=atol
    
    # Test integration of constant function
    @test ∫(0.0, 1.0, x -> 2.0, gausslegendre) ≈ 2.0 atol=atol
    
    # Test integration of sine function
    @test ∫(0.0, π, sin, gausslegendre) ≈ 2.0 atol=atol
    
    # Test integration of exponential function
    @test ∫(0.0, 1.0, exp, gausslegendre) ≈ exp(1.0) - 1.0 atol=atol
end
