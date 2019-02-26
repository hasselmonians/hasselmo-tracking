## Attempt to track a repeating trajectory using the Hille series

# preamble
using LinearAlgebra

# initial vector
x₀ = [0, 1, 2, 3, 0, 1]

# repeat the vector n times
x = repeat(x₀, 3)

# define the forward difference
function forward_diff(x)
    # compute the forward difference of x at ii
    # up to fourth order differences
    @assert length(x) > 5 "length(x) must be > 5"

    coeffs = [-1 1 0 0 0;
              1 -2 1 0 0;
              -1 3 -3 1 0;
              1 -4 6 -4 1]

    dx = [coeffs[i, :] ⋅ x[end-4:end] for i ∈ 1:4]
    return Tuple(dx)
end

# define the Hille series trajectory iteration
function evolve(x, t)
    # iterate x using the forward difference t steps
    y = zeros(length(x) + t)
    dx = forward_diff(x)

    function hille(x, dx, t)
        value = x + t * dx[1] + (1/2)t^2*dx[2] + (1/6)t^3*dx[3] + (1/24)t^4*dx[4]
    end

    for i ∈ 1:t
        y[length(x)+i] = hille(x[end], dx, i)
    end

    return y
end
