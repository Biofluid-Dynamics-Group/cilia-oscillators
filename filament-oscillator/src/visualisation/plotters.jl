using DifferentialEquations
using Plots

"""
    animate_cilia(solution::ODESolution)

Animates the cilia movement given the solution `solution` of the ODE.
"""
function animate_cilia(solution::ODESolution)
    Ψ_array = stack(solution.u, dims=1)[begin:size(solution.u)[1]÷100 + 1:end, :]
    p = plot(
        xlim=(-L*1.1, L*1.1), ylim=(0.0, L*1.1), title="System movement",
        xaxis="x [μm]", yaxis="z [μm]", legend=false
    )
    xlims!(p, (-L*1.1, (M - 1)*d + L*1.1))
    color_scheme = palette(:twilight, size(Ψ_array)[1], rev=true)
    @gif for k = 1:size(Ψ_array)[1]
        x_vector = x(Ψ_array[k, :])
        for j=1:M
            positions = zeros(N+1, 3)
            positions[1, :] = x₀[j]
                for i=1:N
                positions[i+1, :] += x_vector[i + (j - 1)N, :]
            end
            plot!(p, positions[:, 1], positions[:, 3], color=color_scheme[k])
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
