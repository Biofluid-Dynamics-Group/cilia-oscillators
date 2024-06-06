using LinearAlgebra


function oseen_burgers_tensor(x::Vector, y::Vector, μ::Real)
    r = x - y
    return 1/(8*π*μ)*(I/norm(r) + r*r'/norm(r)^3)
end

function potential_dipole(x::Vector, y::Vector, μ::Real)
    r = x - y
    𝛿 = [1 0 0; 0 1 0; 0 0 -1]
    return 1/(8*π*μ)*(I/norm(r)^3 - 3*r*r'/norm(r)^5)*𝛿
end

function stokeslet_doublet(x::Vector, y::Vector, μ::Real)
    r = x - y
    𝛿 = [1 0 0; 0 1 0; 0 0 -1]
    G_D = potential_dipole(x, y, μ)
    return 1/(8*π*μ)*(x[1]*G_D + ((x*[1 0 0] - [1; 0; 0]*x')/(norm(r)^3))*𝛿)
end

function blake_tensor(x::Vector, y::Vector, μ::Real)
    Y = [1, 1, -1].*y
    h = y[3]
    G = oseen_burgers_tensor(x, y, μ)
    Ĝ = oseen_burgers_tensor(x, Y, μ)
    G_D = potential_dipole(x, Y, μ)
    G_SD = stokeslet_doublet(x, Y, μ)
    return G - Ĝ + 2*h^2*G_D - 2*h*G_SD
end

function rotne_prager_blake_tensor(x::Vector, y::Vector, μ::Real, a::Real)
    if norm(x - y) < eps()
        return rotne_prager_blake_self_mobility(x, μ, a)
    else
        return rotne_prager_blake_cross_mobility(x, y, μ, a)
    end
end

function rotne_prager_blake_self_mobility(x::Vector, μ::Real, a::Real)
    X = [1, 1, -1].*x
    v = -3/16*(a/x[3] - (a/x[3])^3 + 1/3*(a/x[3])^5)
    δμ = 1/(6.0*π*μ*a)*[v 0 0; 0 v 0; 0 0 2*v]
    return 1/(6.0*π*μ*a)*I - rotne_prager_tensor(x, X, μ, a) + δμ
end

function rotne_prager_blake_cross_mobility(x::Vector, y::Vector, μ::Real, a::Real)
    Y = [1, 1, -1].*y
    z₁ = x[3]  # Adopting notation from Gauger, Downton & Stark (2009)
    z₂ = y[3]
    s = norm(x - Y)
    R = x - Y  # x, y, z indices for R are 1, 2, 3 respectively
    δμ = zeros(3, 3)
    for α = 1:2
        if α == 1
            β = 2
        else
            β = 1
        end
        δμ[α, α] += -z₁*z₂*(1/s^3 - 3*R[α]^2/s^5)
        δμ[α, α] += -a^2/s^7*R[3]^2*(4*R[α]^2 - R[β]^2 - R[3]^2)
        δμ[α, α] += -a^4/(3*s^9)*(
            4*R[α]^4 - R[β]^4 + 4*R[3]^4 + 3*R[α]^2*R[β]^2 + 3*R[β]^2*R[3]^2 -
            27*R[α]^2*R[3]^2
        )

        δμ[α, β] += 3*z₁*z₂*R[α]*R[β]/s^5
        δμ[α, β] += -5*a^2/s^7*R[α]*R[β]*R[3]^2
        δμ[α, β] += -5*a^4/(3*s^9)*(R[α]^2 + R[β]^2 - 6*R[3]^2)*R[α]*R[β]

        δμ[α, 3] += R[α]*(z₁^2/s^3 - 3*z₁*z₂*R[3]/s^5)
        δμ[α, 3] += R[α]*a^2*(
            1/3*1/s^3 + 5/s^7*R[3]^3 - 1/s^5*(R[3] + 2*(z₁^2 + z₁*z₂))
        )
        δμ[α, 3] += 5*a^4/(3*s^9)*R[α]*(3*R[3]*(R[α]^2 + R[β]^2) - 4*R[3]^3)

        δμ[3, α] += R[α]*(z₁^2/s^3 + 3*z₁*z₂*R[3]/s^5)
        δμ[3, α] += R[α]*a^2*(
            1/3*1/s^3 - 5/s^7*R[3]^3 + 1/s^5*(R[3] - 2*(z₁^2 + z₁*z₂))
        )
        δμ[3, α] += -5*a^4/(3*s^9)*R[α]*(3*R[3]*(R[α]^2 + R[β]^2) - 4*R[3]^3)
    end
    δμ[3, 3] += z₁*z₂*(1/s^3 - 3*R[3]^2/s^5)
    δμ[3, 3] += -a^2/s^7*R[3]^2*(3*(R[1]^2 + R[2]^2) - 2*R[3]^2)
    δμ[3, 3] += -a^4/(3*s^9)*(
        3*(R[1]^4 + R[2]^4) + 6*R[1]^2*R[2]^2 - 24*R[3]^2*(R[1]^2 + R[2]^2) + 8*R[3]^4
    )
    return rotne_prager_tensor(x, y, μ, a) - rotne_prager_tensor(x, Y, μ, a) + δμ
end

function rotne_prager_tensor(x::Vector, y::Vector, μ::Real, a::Real)
    r = x - y
    if norm(r) < eps()
        return 1/(6*π*a*μ)*I
    else
        r_hat = r./norm(r)
        return 1/(6*π*a*μ)*(
            3/4*a/norm(r)*(I + r_hat*r_hat') + 1/2*(a/norm(r))^3*(I - 3*r_hat*r_hat')
        )
    end
end

function RPY_tensor(x::Vector, y::Vector, μ::Real, a::Real)
    r = x - y
    r_norm = norm(r)
    r_hat = r./r_norm
    if r_norm < eps()
        return 1/(6.0*π*μ*a)*I(3)
    else
        return 1/(8.0*π*μ*a)*(
            (1 + 2*a^2/(3*r_norm^2))*I(3) + (1 - 2*a^2/r_norm^2)*r_hat*r_hat'
        )
    end
end
