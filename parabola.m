% based on Mike's "original" code
% seeks to implement the matrix transformation method of following a trajectory

% this code expects only second-order accuracy
% otherwise the matrices would have to be bigger
% e.g. 4 x 4 matrices for third-order accuracy

qqq

%% Instantiate important variables

% time resolution
dt      = 0.01;
t_end   = 10;
time    = dt:dt:t_end;
nSteps  = length(time);
mindex  = 6;

% generate the trajectory
x       = 10*time;
y       = -40 + 100*time - 10*time.^2;

% generate the corner matrix
c       = ones(3);

for ii = 1:mindex
  c0    = c;
  

%%%%%%%

% trajectory-tracking matrices
corner  = ones(3);
W       = NaN(2, nSteps, 3, 3);
  % the first dimension is velocity or acceleration
  % the second dimension is the time index
  % the third and fourth dimensions produce the 2-D matrix
xmat    = NaN(mindex, 3);
ymat    = NaN(mindex, 3);
xymat   = ones(2, 2, 3, 3);
  % the first dimension is x or y
  % the second dimension is which time index it is
  % the third and fourth dimensions produce the 2-D matrix
Wx      = NaN(mindex, 3, 3);
Wy      = NaN(mindex, 3, 3);
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

for ii = 1:mindex
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
    W(2, ii, :, :)  = squeeze(W(1, ii, :, :)) * inv(squeeze(W(1, ii-1, :, :)));

    % store the position, velocity, and acceleration in matrices
    xmat(ii, :)     = [x(ii), W(1, ii, 1, 3), W(2, ii, 1, 3)];
    ymat(ii, :)     = [y(ii), W(1, ii, 2, 3), W(2, ii, 2, 3)];
  end

  if ii > 3
    % create sequential position, velocity, acceleration arrays
    % restructure xmat and ymat into a 3x3 ma                                                                                                                    trix
    xymat(1, 1, 1:3, 1:2) = xmat(ii-3:ii-1, 1:2);
    xymat(1, 2, 1:3, 1:2) = xmat(ii-2:ii, 1:2);
    xymat(2, 1, 1:3, 1:2) = ymat(ii-3:ii-1, 1:2);
    xymat(2, 2, 1:3, 1:2) = ymat(ii-2:ii, 1:2);

    % compute the inverse prediction matrix
    Wx(ii, :, :)    = squeeze(xymat(1, 2, :, :)) * inv(squeeze(xymat(1, 1, :, :)));
    Wy(ii, :, :)    = squeeze(xymat(2, 2, :, :)) * inv(squeeze(xymat(2, 1, :, :)));
  end
end

%% Predict the Trajectory

for ii = 7:nSteps
  xhat(ii, :)       = squeeze(Wx(mindex, :, :))^(ii-mindex) * transpose([xmat(mindex, 1:2), 1]);
  yhat(ii, :)       = squeeze(Wy(mindex, :, :))^(ii-mindex) * transpose([ymat(mindex, 1:2), 1]);
end

%% Plot the Trajectory

figure;
plot(x, y, 'ok')
plot(xhat(:, 1), yhat(:, 1), 'or')
