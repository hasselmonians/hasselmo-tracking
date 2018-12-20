%% Hille Series Trajectory Tracing
% Alec Hoyland
% 2018-12-20 10:59

pHeader;
tic

%% Introduction
% The Hille series is equivalent to a discretized Taylor series under the limit
% $$ \lim_{\Delta t \rightarrow 0} \sum_{n=0}^{\infty} \frac{t^n}{n! (\Delta t)^n} \mathrm{D}^n f(a) = f(a + t) $$
% for $t > 0$ and $\mathrm{D}^n$ is the finite difference operator of order $n$.
%
% For a discrete time step $\Delta t$, the trajectory $f$ can be predicted at future times $a + t$.
% The number of historical trajectory points needed depends on the order of the approximation.
% When expanded, this equation yields:
% $$ \left[ 1 + \frac{t}{\Delta t} \mathrm{D}^1 + \frac{t^2}{2 (\Delta t)^2} \mathrm{D}^2 + ... \right] f(a) $$
%
% The finite difference operator combines past trajectory points in proportions according to order.
% The first three terms of the Hille series in matrix form are:
% $$ \pmatrix{1 \frac{t}{\Delta t} \frac{t^2}{2(\Delta t)^2}}
% \pmatrix{0 0 1 \cr -\frac{1}{2} 0 \frac{1}{2} \cr 1 -1 1 \cr}
% \pmatrix{f(a) \cr f(a-\Delta t) \cr f(a - 2\Delta t)} $$

%% Version Info
pFooter;

t = toc;

%%
% This document was built in:
disp(strcat(oval(t,3),' seconds.'))
