using Test

include("../src/models/cilia.jl")

fluid = FluidParameters(μ)
@testset "Identical forcings" begin
    for M_test = 1:5
        for ψ₁ = 0.0:0.1:2π
            h = L/(N - 1)
            params = CiliaParameters(M_test, N, φ, d, a, κ_b)
            x₀ = [[(j - 1)*d, 0.0, 0.0] for j=1:M]
            s = [(i)*h for i=1:N]
            system = CiliaSystem(
                params,
                x₀,
                s,
                I(3),
                repeat([ψ₁, 0.0], outer=[M_test])
            )

            Q = Q_ref(system, fluid)
            K_h_matrix = convert(Matrix{Float64}, Kₕ(1, system))
            q = K_h_matrix'*(Mₕ(1, system, fluid)\K_h_matrix*[ω, 0.0])
            for j = 2:M_test
                @test Q[2j - 1:2j] ≈ q atol=1e-10
            end
        end
    end
end
