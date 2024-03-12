using Plots

include("Rotors.jl")
include("DefaultParameters.jl")


λ_array = [0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 100].*λ
# λ_array = [100*λ]
num_steps = 100000

for λi in λ_array
    sol = simulate_rotor(
        (λ=λi, η=η, f_drive=f_drive, r_0=r_0, a=a, x_0=x_0, gamma_0=gamma_0),
        final_time, num_steps
    )
    x_vals = sol.u
    t_vals = sol.t
    last_periodish_index = findall(t_vals .> final_time-2*T)[1] - 1
    position_periodish = x_vals[last_periodish_index:end]
    time_periodish = t_vals[last_periodish_index:end]

    r_vals = zeros(length(time_periodish))
    for i in eachindex(r_vals)
        r_vals[i] = r(position_periodish[i], x_0)
    end

    position_periodish = stack(position_periodish)

    r_vals = (r_vals .- r_0) ./ r_0

    display(plot!(time_periodish, r_vals))
    # display(plot!(position_periodish[1, begin:end], position_periodish[3, begin:end]))
    # display(plot!(sol, vars=(1, 3), aspect_ratio=:equal))
end
