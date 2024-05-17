using LinearAlgebra


A_1 = [
    -0.327 0.393 -0.097 0.079;
    0.3935 -1.516 0.032 -0.302;
    0.101 0.716 -0.118 0.142
]
A_2 = [
    0.9475 -0.018 0.158 0.01;
    -0.276 -0.126 -0.341 0.035;
    0.048 0.263 0.186 -0.067
]
B_1 = [
    0.0 0.284 0.006 -0.059;
    0.0 1.045 0.317 0.226;
    0.0 -1.017 -0.276 -0.196
]
B_2 = [
    0.0 0.192 -0.05 0.012;
    0.0 -0.499 0.423 0.138;
    0.0 0.339 -0.327 -0.114
]

M = 3
N = 3

function ξ(s::Real, ϕ::Real)
    return [dot(s.^collect(1:M), A_1*cos.(collect(0:N)*ϕ) + B_1*sin.(collect(0:N)*ϕ)),
    dot(s.^collect(1:M), A_2*cos.(collect(0:N)*ϕ) + B_2*sin.(collect(0:N)*ϕ))]
end

function dξ_dϕ(s::Real, ϕ::Real)
    return [dot(s.^collect(1:M), -A_1*sin.(collect(0:N)*ϕ) + B_1*cos.(collect(0:N)*ϕ)),
    dot(s.^collect(1:M), -A_2*sin.(collect(0:N)*ϕ) + B_2*cos.(collect(0:N)*ϕ))]
end    
