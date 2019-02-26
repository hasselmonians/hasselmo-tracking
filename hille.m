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
% $$ \pmatrix{1 \quad \frac{t}{\Delta t} \quad \frac{t^2}{2(\Delta t)^2}}
% \pmatrix{0 \quad 0 \quad 1 \cr -\frac{1}{2} \quad 0 \quad \frac{1}{2} \cr 1 \quad -1 \quad 1 \cr}
% \pmatrix{f(a) \cr f(a-\Delta t) \cr f(a - 2\Delta t)} $$

%% An Example
% The trajectory we will consider is a simple polynomial.

% initialize the time step
dt    = 0.1; % s

% the perceiver has seen 1 second of data
a     = 0:dt:1;

% define the feature matrix
fcn   = @(x) transpose(3*x - x.^2);
f     = fcn(a);

% compute the Hille coefficients
D     = getFiniteDifferenceCoeffs(6);
% compute the kernel
k     = D * f(end:-1:end-6);

return

% set up the real trajectory
t     = dt:dt:1;
traj0 = fcn([a, a(end) + t]); % the real trajectory

% iterate the model forwards for another second
traj  = [fcn(a); NaN(length(t), 1)]; % the predicted trajectory
T     = zeros(1, 7); % stores the prefactors
for ii = 1:length(t)
  for ww = 1:7
    T(ww)   = t(ii)^(ww-1) * prefactor(ww-1, dt);
  end
  qq        = ii + length(a);
  traj(qq)  = T * k;
end

figure('OuterPosition',[0 0 1200 1200],'PaperUnits','points','PaperSize',[1200 1200]); hold on
plot([a, a(end) + t], traj0, 'k');
plot([a, a(end) + t], traj, 'r');
xlabel('time (s)')
ylabel('f(t)')
legend({'trajectory', 'prediction'}, 'Location', 'best')
title('Predicting a sinusoidal trajectory')

figlib.pretty()

if being_published
  snapnow
  delete(gcf)
end


%% Version Info
pFooter;

timeDocumentWasBuiltIn = toc;

%%
% This document was built in:
disp(strcat(oval(timeDocumentWasBuiltIn,3),' seconds.'))

function y = prefactor(n, dt)
  % accepts the term order, the time step, and the forward time (in units of time)
  % and computes the prefactor to the nth term in the Hille series
  y = 1 / (factorial(n) * dt^n);
end % prefactor

function coeff = getFiniteDifferenceCoeffs(n)
  % n is the order of derivative you want

  if n == 0
    coeff = [1];
  elseif n == 1
    coeff = [1, 0, 0; 1/2, 0, -1/2];
  elseif n == 2
    coeff = [1, 0, 0; 1/2, 0, -1/2; -1, 2, -1];
  elseif n == 3
    coeff = [1, 0, 0, 0, 0;, 1/2, 0, -1/2, 0, 0; -1, 2, -1, 0, 0; 1/2, -1, 0, 1, -1/2];
  elseif n == 4
    coeff = [1, 0, 0, 0, 0; ...
             1/2, 0, -1/2, 0, 0; ...
             -1, 2, -1, 0, 0; ...
             1/2, -1, 0, 1, -1/2; ...
             1, -4, 6, -4, 1];
  elseif n == 5
    coeff = [1, 0, 0, 0, 0, 0, 0; ...
             1/2, 0, -1/2, 0, 0, 0, 0; ...
             -1, 2, -1, 0, 0, 0, 0; ...
             1/2, -1, 0, 1, -1/2, 0, 0; ...
             1, -4, 6, -4, 1, 0, 0; ...
             1/2, -2, 5/2, 0, -5/2, 2, -1/2];
  elseif n == 6
     coeff = [1, 0, 0, 0, 0, 0, 0; ...
              1/2, 0, -1/2, 0, 0, 0, 0; ...
              -1, 2, -1, 0, 0, 0, 0; ...
              1/2, -1, 0, 1, -1/2, 0, 0; ...
              1, -4, 6, -4, 1, 0, 0; ...
              1/2, -2, 5/2, 0, -5/2, 2, -1/2; ...
              1, -6, 15, -20, 15, -6, 1];
  else
    coeff = zeros(n);
    coeff(1:7, 1:7) = getFiniteDifferenceCoeffs(6);
  end

  % coeff = -coeff;

end % getFiniteDifferenceCoeffs
