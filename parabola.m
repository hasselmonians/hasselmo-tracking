% based on Mike's "original" code
% seeks to implement the matrix transformation method of following a trajectory

% this code expects only second-order accuracy
% otherwise the matrices would have to be bigger
% e.g. 4 x 4 matrices for third-order accuracy

qqq

%% Instantiate important variables

% time resolution
dt      = 0.1;
t_end   = 10;
time    = dt:dt:t_end;
nSteps  = length(time);

% trajectory-tracking matrices
corner  = ones(3);
W       = NaN(2, nSteps, 3, 3);
  % the first dimension is velocity or acceleration
  % the second dimension is the time index
  % the third and fourth dimensions produce the 2-D matrix
xmat    = NaN(6, 3);
ymat    = NaN(6, 3);
xymat   = ones(2, 2, 3, 3);
  % the first dimension is x or y
  % the second dimension is which time index it is
  % the third and fourth dimensions produce the 2-D matrix
Wx      = NaN(6, 3, 3);
Wy      = NaN(6, 3, 3);
xhat    = NaN(nSteps, 3);
yhat    = NaN(nSteps, 3);

% generate the trajectory
x       = 10*time;
y       = -40 + 100*time - 10*time.^2;

%% Generate prediction matrices
% update corner that reflects three features
% based on position, velocity

% time-evolution matrix
tstepmat          = zeros(3);
tstepmat(1, 1)    = dt;
tstepmat(2, 3)    = dt;

for ii = 1:6
  % save the old corner
  oldcorner       = corner;

  % build the basic corner at t(ii+1)
  corner(1, :)    = x(ii+1);
  corner(2, :)    = y(ii+1);
  corner(3, :)    = 1;

  % add the time-evolution matrix
  corner          = corner + tstepmat;

  if ii > 1
    % compute the transformation matrix for velocity
    W(1, ii, :, :)  = corner * inv(oldcorner);
  end

  if ii > 2
    % compute the transformation matrix for acceleration
    % index by time and reduce to 2 dimensions
    WW              = squeeze(W(1, ii, :, :));
    WW0             = squeeze(W(1, ii-1, :, :));

    % compute the acceleration matrix and store
    W(2, ii, :, :)  = WW * inv(WW0);

    % store the position, velocity, and acceleration in matrices
    xmat(ii, :)     = [x(ii) W(1, ii, 1, 3) W(2, ii, 1, 3)];
    ymat(ii, :)     = [y(ii) W(1, ii, 2, 3) W(2, ii, 2, 3)];
  end

  if ii > 3
    % create sequential position, velocity, acceleration arrays
    % restructure xmat and ymat into a 3x3 matrix
    xymat(1, 1, 1:3, 1:2) = xmat(ii-3:ii-1, 1:2);
    xymat(1, 2, 1:3, 1:2) = xmat(ii-2:ii, 1:2);
    xymat(2, 1, 1:3, 1:2) = ymat(ii-3:ii-1, 1:2);
    xymat(2, 2, 1:3, 1:2) = ymat(ii-2:ii, 1:2);

    % compute the inverse prediction matrix
    % for x-dimension
    WW0             = squeeze(xymat(1, 1, :, :));
    WW              = squeeze(xymat(1, 2, :, :));
    Wx(ii, :, :)    = WW * inv(WW0);

    % for y-dimension
    WW0             = squeeze(xymat(2, 1, :, :));
    WW              = squeeze(xymat(2, 2, :, :));
    Wy(ii, :, :)    = WW * inv(WW0);
  end
end

%% Predict the Trajectory

for ii = 7:nSteps
  xhat(ii, :)       = squeeze(Wx(6, :, :))^(ii-6) * transpose([xmat(6, 1:2), 1]);
  yhat(ii, :)       = squeeze(Wy(6, :, :))^(ii-6) * transpose([ymat(6, 1:2), 1]);
end

%% Plot the Trajectory

figure;
plot(x, y, 'ok')
plot(xhat(:, 1), yhat(:, 1), 'or')
