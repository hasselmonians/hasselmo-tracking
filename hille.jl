## The Hille Series

# time resolution
dt          = 0.1; # seconds
t_end       = 20; # seconds
t           = 0:dt:t_end;

# feature dynamics
# the trajectory is going to be a quadratic with a 1st order term
function trajectory(t, params)
    # t is the time at which the trajectory should be computed
    # params is the set of three parameters (x_0, v_0, a_0)
    return params[1] + params[2] * t + 0.5 * params[3] * t^2;
end

# define the coefficients for a 2nd order fit
coeff   = [1 0 0; 1/2 0 -1/2; -1 2 -1;];

# define the prefactor
function prefactor!(P::Array{T, 2}, t::T, dt::T) where T <: Real
    # assume second order
    P[1]    = 1;
    P[2]    = t / dt;
    P[3]    = t^2 / dt^2 / 2.0;
    return P
end

# construct the real trajectory
x       = zeros(length(t));
x       = [trajectory(i, (1, 100, -10)) for i in t]

# compute the predicted trajectory
function computeTrajectory(x, t, dt)
    # assuming a second-order fit
    if length(x) != length(t)
        error("x and t are not the same length")
    end
    # set up the output vector with known values
    x′      = zeros(length(x), 1);
    x′[1:3] = x[1:3];
    # set up the prefactor vector
    P       = zeros(Float64, 1, 3);
    temp    = zeros(Float64, 1, 1);
    for ii in 4:length(t)
        # compute the prefactor
        prefactor!(P, t[ii], dt);
        # compute the predicted trajectory at time t[ii]
        temp = P * coeff * x′[ii-3:ii-1, :];
        x′[ii] = reshape(temp, 1)[1];
    end
    return x′
end # computeTrajectory

# plot to check the trajectory
using Plots; plotlyjs()
plot(t, x)
x′ = computeTrajectory(x, t, dt)
plot!(t, x′)
