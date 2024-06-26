using Base.Threads
using LinearAlgebra

include("../beats/tangent_angle.jl")
include("../utils/greens_functions.jl")
include("../utils/algebra.jl")


struct FluidParameters{T<:Real}
    μ::T  # Viscosity
end

struct SimulationParameters{FloatType<:Real, IntType<:Integer}
    M::IntType      # Number of cilia
    N::IntType      # Number of discretisation spheres per cilium
    φ::FloatType    # Angle of cilia beat plane
    d::FloatType    # Distance between cilia
    a::FloatType    # Radius of discretisation spheres
    κ_b::FloatType  # Bending stiffness
end

struct CiliaSystem{T<:Real}
    beat_params::BeatParameters       # Beat parameters
    sim_params::SimulationParameters  # Simulation parameters
    x₀::Vector{Vector{T}}             # Basal positions of cilia
    s::Vector{T}                      # Arclength position of discretisation spheres
    A::Matrix{T}                      # Beat plane rotation matrix
    Ψ::Vector{T}                      # System phase
    u::Function                       # Background flow
end

"""
    x(system::CiliaSystem)

Returns the stacked positions of the discretised cilia. It is an NM x 3 size matrix
with the positions of the `j`th cilium at the `i`th discretisation point given by
`x[i + (j - 1)N, :]` where each cilium is discretised into N spheres.
"""
function x(system::CiliaSystem)
    num_positions = system.sim_params.M*system.sim_params.N
    positions = zeros(3, num_positions)  # Transposed for efficient memory access
    @threads for idx=1:num_positions
        j = ceil(Int, idx / system.sim_params.N)
        i = idx - (j - 1)*system.sim_params.N
        positions[:, idx] .= system.A*ξ(
            system.s[i], system.Ψ[2j - 1:2j], system.beat_params
        ) + system.x₀[j]
    end
    return positions'
end

"""
    Kₕ(j::Int, system::CiliaSystem)

Returns the matrix that maps the phase velocity of the `j`th cilium to its velocities.
"""
function Kₕ(j::Int, system::CiliaSystem)
    # Each matrix is the stacking of the individual K matrices along the cilium length
    return cat(
        [system.A*K(
            system.s[i], system.Ψ[2j - 1:2j], system.beat_params
        ) for i = 1:system.sim_params.N]..., dims=1
    )
end

"""
    κₕ(system::CiliaSystem)

Returns the block-diagonal matrix that maps the phase velocity to the discretised cilia
system velocity.
"""
function κₕ(system::CiliaSystem)
    # Initialise array to hold the Kₕ matrices for each cilium
    K_h_matrix_array = zeros(system.sim_params.M, 3system.sim_params.N, 2)
    # Populate the array. Cilia are geometrically independent, so threads can be used
    for j=1:system.sim_params.M  # Inefficient looping...
        # Each matrix is the stacking of the individual K matrices along the cilium length
        K_h_matrix_array[j, :, :] .= Kₕ(j, system)
    end
    return block_diagonal(K_h_matrix_array)
end

"""
    Πₕ(system::CiliaSystem, fluid::FluidParameters)

Returns the block-diagonal regularised Stokes mobility matrix mapping forces in the
discretised cilia system to velocities.
"""
function Πₕ(system::CiliaSystem, fluid::FluidParameters)
    num_positions = system.sim_params.N*system.sim_params.M
    mobility = zeros(3*num_positions, 3*num_positions)
    x_vector = x(system)
    # Forcings depend on positions which are independent, so threads can be used
    for j = 1:num_positions
        for i = 1:num_positions
            α = 1 + 3*(i - 1)
            β = 1 + 3*(j - 1)
            tensor = RPY_tensor(
                x_vector[i, :], x_vector[j, :], fluid.μ, system.sim_params.a
            )
            mobility[α:(α + 2), β:(β + 2)] .= tensor
        end
    end
    return mobility
end

"""
    Mₕ(j::Int, system::CiliaSystem, fluid::FluidParameters)

Returns the regularised Stokes mobility matrix mapping forces for the `j`th cilium in the
system to its velocities. Used to calculate independent forcings.
"""
function Mₕ(k::Int, system::CiliaSystem, fluid::FluidParameters)
    mobility = zeros(3*system.sim_params.N, 3*system.sim_params.N)
    ψ = system.Ψ[2k - 1:2k]  # Specific cilium phase
    for j = 1:system.sim_params.N
        for i = 1:system.sim_params.N
            α = 1 + 3*(i - 1)
            β = 1 + 3*(j - 1)
            mobility[α:(α + 2), β:(β + 2)] .= RPY_tensor(
                    system.A*ξ(system.s[i], ψ, system.beat_params),
                    system.A*ξ(system.s[j], ψ, system.beat_params),
                    fluid.μ, system.sim_params.a
                )
        end
    end
    return mobility
end

"""
    Q_ref(system::CiliaSystem, fluid::FluidParameters)

Returns the generalised force vector that, in absence of other cilia, induces a
reference beat in each filament.
"""
function Q_ref(system::CiliaSystem, fluid::FluidParameters)
    forcings = zeros(2system.sim_params.M)
    # For each individual cilium, calculate the forcing that induces a reference beat
    # independently
    for j=1:system.sim_params.M
        K_h_matrix = convert(Matrix{Float64}, Kₕ(j, system))
        M_h_matrix = Mₕ(j, system, fluid)
        forcings[2j - 1:2j] .= K_h_matrix'*inv(M_h_matrix)*K_h_matrix*[
            2π/system.beat_params.T, 0.0
        ]
    end
    return forcings
end
