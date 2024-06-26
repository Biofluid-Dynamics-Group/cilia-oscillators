using DifferentialEquations
using Plots

include("../beats/tangent_angle.jl")
include("../models/cilia.jl")	


"""
    animate_cilia(solution::ODESolution)

Animates the cilia movement given the solution `solution` of the ODE.
"""
function animate_cilia(solution::ODESolution, system::CiliaSystem)
    Ψ_array = stack(solution.u, dims=1)[begin:size(solution.u)[1]÷200 + 1:end, :]
    p = plot(
        ylim=(0.0, L*1.1), title="System movement",
        xaxis="x [μm]", yaxis="z [μm]", legend=false
    )
    xlims!(p, (-L*1.1, (M - 1)*d + L*1.1))
    color_scheme = palette(:twilight, size(Ψ_array)[1], rev=true)
    @gif for k = 1:size(Ψ_array)[1]
        system.Ψ .= Ψ_array[k, :]
        x_vector = x(system)
        for j=1:system.sim_params.M
            positions = zeros(N+1, 3)
            positions[1, :] = system.x₀[j]
                for i=1:system.sim_params.N
                positions[i+1, :] += x_vector[i + (j - 1)N, :]
            end
            plot(p, positions[:, 1], positions[:, 3], color=color_scheme[k])
        end
    end
    display(p)
    return nothing
end

"""
    plot_system_trajectory(solution::ODESolution)

Plots the phase evolution of the system given the solution `solution` of the ODE.
"""
function plot_system_trajectory(solution::ODESolution)
    Ψ_array = stack(solution.u, dims=1)
    t = solution.t
    p = plot(title="Phase evolution", xaxis="t [ms]")
    for i=1:size(Ψ_array)[2]
        plot!(t, Ψ_array[:, i], label="ψ_$i")
    end
    display(p)
    return nothing
end

function plot_at_phase(Ψ::Vector, system::CiliaSystem)
    system.Ψ .= Ψ
    x_vector = x(system)
    p = plot()
    for j=1:system.params.M
        positions = zeros(N+1, 3)
        positions[1, :] = system.x₀[j]
            for i=1:system.params.N
                positions[i+1, :] += x_vector[i + (j - 1)N, :]
            end
        plot!(p, positions[:, 1], positions[:, 3])
    end
    display(p) 
end


# looking at Ψ when the system fails, is it consistent?

function animate_cilia(data::DataFrame, system::CiliaSystem)
    xlimit = system.beat_params.L*1.1*(system.sim_params.M)
    p = plot(
        ylim=(0.0, system.beat_params.L*1.1), xlim=(-system.beat_params.L*1.1, xlimit),title="System movement",
        xaxis="x [μm]", yaxis="z [μm]", legend=false
        )
    color_scheme = palette(:twilight, size(data)[1], rev=true)
    @gif for k = 1:size(data)[1]
        system.Ψ .= data[k, :Ψ]
        x_vector = x(system)
        for j=1:system.sim_params.M
            positions = zeros(N+1, 3)
            positions[1, :] = system.x₀[j]
            for i=1:system.sim_params.N
                positions[i+1, :] += x_vector[i + (j - 1)N, :]
            end
            plot(p, positions[:, 1], positions[:, 3], color=color_scheme[k])
        end
    end
end
