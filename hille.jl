## The Hille Series
# predict a trajectory by using the Hille series

# preamble
using Calculus, Plots; plotlyjs();

# instantiate time parameters
Δt      = 0.01
t_end   = 3
t       = Δt:Δt:t_end

# generate the time-series
f(t::Real)  = -40 + 100t - 10t^2
ḟ(t::Real)  = Calculus.derivative(f, t)
f̈(t::Real)  = Calculus.second_derivative(f, t)

x           = f.(t)
v           = ḟ.(t)
a           = f̈.(t)

# compute the approximation using the Taylor series
xT          = zeros(length(t))
Taylor(t, x, v, a) = x + v*t + 0.5*a*t^2
(x0, v0, a0) = (x[3], v[3], a[3])
for ii in 4:length(t)
    xT[ii]  = Taylor(t[ii], x0, v0, a0)
end

## Plot the results
plot(t, x)
plot!(t, xT)
