%% Hille Series Trajectory Tracing
% Alec Hoyland
% 2018-12-20 10:59

pHeader;
tic

%% Introduction
% The Hille series is equivalent to a discretized Taylor series under the limit
% $$ \lim_{\Delta t \rightarrow 0} \sum_{n=0}^{\infty} \frac{t^n}{n! \Delta t^n} \mathrm{D}^n f(a) = f(a + t)$$ for $t > 0$.

%% Version Info
pFooter;

t = toc;

%%
% This document was built in:
disp(strcat(oval(t,3),' seconds.'))
