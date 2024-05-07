using LinearAlgebra
using Base.Threads


"""
    block_diagonal(matrices::Array{Matrix})

Given an array of equally-sized matrices, returns a block-diagonal matrix with the
matrices of the array on the diagonal.
"""
function block_diagonal(matrices::Vector)
    example = matrices[1]
    rows, columns = size(example)
    n = length(matrices)
    result = zeros(rows*n, columns*n)
    @threads for i = 1:n
        result[1 + (i - 1)*rows:i*rows, 1 + (i - 1)*columns:i*columns] .= matrices[i]
    end
    return result
end
function block_diagonal(matrices::Array{Float64, 3})
    example = matrices[1, :, :]
    rows, columns = size(example)
    n = size(matrices)[1]
    result = zeros(rows*n, columns*n)
    @threads for i = 1:n
        result[1 + (i - 1)*rows:i*rows, 1 + (i - 1)*columns:i*columns] .= matrices[i, :, :]
    end
    return result
end
