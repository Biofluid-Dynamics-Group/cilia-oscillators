using Base.Threads
using LinearAlgebra

include("../beats/tangent_angle.jl")
include("../utils/greens_functions.jl")


struct FluidParameters
    μ::Real  # Viscosity
end

struct CiliaParameters
    M::Int  # Number of cilia
    N::Int  # Number of discretisation spheres per cilium
    φ::Real  # Angle of cilia beat plane
    d::Real  # Distance between cilia
    a::Real  # Radius of discretisation spheres
    κ_b::Real  # Bending stiffness
end

struct CiliaSystem
    params::CiliaParameters  # System parameters
    x₀::Vector{Vector{Real}}  # Basal positions of cilia
    s::Vector{Real}  # Arclength position of discretisation spheres
    A::Matrix{Real}  # Beat plane rotation matrix
    Ψ::Vector{Float64}  # System phase
end

"""
    x(system::CiliaSystem, params::CiliaParameters)

Returns the stacked positions of the discretised cilia. It is an NM x 3 size matrix
with the positions of the `j`th cilium at the `i`th discretisation point given by
`x[i + (j - 1)N, :]` where each cilium is discretised into N spheres.
"""
function x(system::CiliaSystem)
    num_positions = system.params.M*system.params.N
    positions = zeros(num_positions, 3)
    @threads for idx=1:num_positions
        j = ceil(Int, idx / system.params.N)
        i = idx - (j - 1)*system.params.N
        positions[idx, :] .= system.A*ξ(system.s[i], system.Ψ[2j - 1:2j]) + system.x₀[j]
    end
    return positions
end

"""
    Kₕ(j::Int, system::CiliaSystem)

Returns the matrix that maps the phase velocity of the `j`th cilium to its velocities.
"""
function Kₕ(j::Int, system::CiliaSystem)
    # Each matrix is the stacking of the individual K matrices along the cilium length
    return cat(
        [system.A*K(system.s[i], system.Ψ[2j - 1:2j]) for i = 1:params.N]..., dims=1
    )
end

"""
    κₕ(system::CiliaSystem)

Returns the block-diagonal matrix that maps the phase velocity to the discretised cilia
system velocity.
"""
function κₕ(system::CiliaSystem)
    # Initialise array to hold the Kₕ matrices for each cilium
    K_h_matrix_array = zeros(system.params.M, 3system.params.N, 2)
    # Populate the array. Cilia are geometrically independent, so threads can be used
    @threads for j=1:system.M
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
    num_positions = system.params.N*system.params.M
    mobility = zeros(3*num_positions, 3*num_positions)
    x_vector = x(system)
    # Forcings depend on positions which are independent, so threads can be used
    @threads for i = 1:num_positions
        @threads for j = 1:num_positions
            α = 1 + 3*(i - 1)
            β = 1 + 3*(j - 1)
            tensor = rotne_prager_blake_tensor(
                x_vector[i, :], x_vector[j, :], fluid.μ, system.a
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
function Mₕ(j::Int, system::CiliaSystem, fluid::FluidParameters)
    mobility = zeros(3*system.params.N, 3*system.params.N)
    ψ = system.Ψ[2j - 1:2j]  # Specific cilium phase
    for i = 1:params.N
        for j = 1:params.N
            α = 1 + 3*(i - 1)
            β = 1 + 3*(j - 1)
            mobility[α:(α + 2), β:(β + 2)] .= rotne_prager_blake_tensor(
                    ξ(system.s[i], ψ), ξ(system.s[j], ψ), fluid.μ, system.params.a
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
    forcings = zeros(2system.params.M)
    # For each individual cilium, calculate the forcing that induces a reference beat
    # independently
    @threads for j=1:system.params.M
        K_h_matrix = Kₕ(j, system)
        forcings[2j - 1:2j] .= K_h_matrix'*(Mₕ(j, system, fluid)\K_h_matrix*[2π, 0.0])
    end
    return forcings
end
