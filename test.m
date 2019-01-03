%% Test the accuracy of 2-D path integration
% comparing the matrix transformation model to a physical model

pHeader;

%% Matrix transformation model
% The "matrix transformation" model describes evolution over time of a 2-D
% trajectory by describing the contributions to position, velocity, and acceleration
% as elements of a matrix. The matrix is $ 3 \times 3$ in order to account for the 0th,
% 1st, and 2nd time derivatives.

% time resolution
tEnd        = 50;
dt          = 1;
nSteps      = tEnd / dt;

% instantiate x and y position, velocity, acceleration vectors
xmat        = zeros(3, nSteps);
ymat        = zeros(3, nSteps);

% transition matrices
Wx          = [1, 1, 0; 0, 1, 1; 0, 0, 1];
Wy          = [1, 1, 1; 0, 1, 1; 0, 0, 1];

disp('The x-dimension transition matrix')
disp(Wx)

disp('The y-dimension transition matrix')
disp(Wy)

% initial conditions
xmat(:, 1)  = transpose([1, 1, 0]);
ymat(:, 1)  = transpose([2, 20, -0.8]);

%% Model description
% We will consider a quadratic trajectory of a particle with a constant acceleration
% in the y-direction and zero acceleration in the x-direction. Both directions have
% positive initial positions and velocities.

% perform simulation using matrix transformation formalism
for ii = 1:nSteps-1
    xmat(:, ii+1)   = Wx * xmat(:, ii);
    ymat(:, ii+1)   = Wy * ymat(:, ii);
end


%% A physical description
% The system depicted here is of a two-dimensional trajectory with constant acceleration.
% As such, the trajectory can be computed as $d = d_0 + v t + \frac{1}{2} a t^2$.

t       = 0:dt:tEnd-dt;
x       = zeros(1, length(t));
y       = zeros(1, length(t));
x       = xmat(1, 1) + xmat(2, 1) * t + (1/2) * xmat(3, 1) * t.^2;
y       = ymat(1, 1) + ymat(2, 1) * t + (1/2) * ymat(3, 1) * t.^2;

% plot the result
figure; hold on
plot(xmat(1, :), ymat(1, :), 'k');
plot(x, y, 'r');

xlabel('x-position')
ylabel('y-position')
title('quadratic trajectory')
legend({'matrix model', 'physical model'})

prettyFig()

if being_published
    snapnow
    delete(gcf)
end

%% A minor adjustment
% Note that if we change the transition matrix in $y$ to

Wy(1, 3)    = 0.5;

disp(Wy)

%%
% then we will see perfect accordance between the theoretical "physical" model
% and the matrix transformation model. This is because while the matrix transformation
% model is an accurate representation of the physical model in matrix form, it must
% rely on knowing the functional form of the trajectory _a priori_, just as the
% physical model must.

pFooter;
