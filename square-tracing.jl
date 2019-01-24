# trace a square that increases in side length with time

# preamble
using Calculus

# feature evolution rule
function feature(t)
    # the eight dimensions correspond to the x- and y-coordinates
    # of the vertices of the square
    x       = zeros(8)
    x[1]    = 0 # top left x
    x[2]    = 0 # top left y
    x[3]    = t # top right x
    x[4]    = 0 # top right y
    x[5]    = 0 # bottom left x
    x[6]    = t # bottom left y
    x[7]    = t # bottom right x
    x[8]    = t # bottom right y
    return x
end

# true feature vector
function feature_evolve(t)
    f   = zeros(length(t), 8)
    for tt in t
        f[tt, :] = feature(tt)
    end
    return f
end

# approximate feature vector
feature_diff(t) = Calculus.derivative(feature, t)
feature_approx(t) = zeros(8) + feature_diff(t)*t

function feature_approx_evolve(t)
    f = zeros(length(t), 8)
    for tt in t
        f[tt, :] = feature_approx(tt)
    end
    return f
end

# instantiate time vector
t       = 1:5

# generate the feature vectors
f       = feature_evolve(t)
f̂       = feature_approx_evolve(t)

# will be true if the approximation is good
isapprox(f, f̂)
