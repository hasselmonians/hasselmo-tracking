# an attempt to represent a change in trajectories using the Laplace transform

using DSP, Plots

## Define the signal
# let the signal be a sinusoid
# the rule is that the sinusoid increases in amplitude by 20% and frequency by 20%
# at the first step, it is a sinusoid: x(t) = A*sin(w*t) where A = w = 1

# set up the time vector
t     = 0:0.01:6π;
timeseries(t, s) = s * sin.(s * t);

## Given the first two patterns in the series
# compute the difference using the Laplace transform

x1 = timeseries(t, 1);
x2 = timeseries(t, 1.2);
Δx = 
