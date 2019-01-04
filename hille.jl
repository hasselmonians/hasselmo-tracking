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

# construct the real trajectory
x       = zeros(length(t));
x       = [trajectory(i, (1, 100, -10)) for i in t]

#
