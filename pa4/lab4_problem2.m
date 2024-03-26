%Problem Setup
%pendulum = Pendulum();
mass = pendulum.mass;
g = pendulum.g;
length = pendulum.length_to_COM;
delta_t = 0.001; % Set the time-step
% t_final = 10.0; % Set the final time.
% Create an evenly spaced grid of times from
% 0 to t_final using the time step delta_t.
t_grid = 0:delta_t:t_final;
% Set N to the number of time steps.
N = numel(t_grid)
% Set the initial condition (theta0 and omega0).
x0 = [pi; 4]; 
% Generate data
[theta_data , omega_data] ...
= pendulum.simulateAndMeasure(t_grid , x0);

disp("2.1.1");
% Pendulum properties
mass = 1.0; % Mass of the pendulum in kg
g = 9.8; % Acceleration due to gravity in m/s^2 (negative indicates direction)
length_to_COM = 1.0; % Length from pivot to center of mass in meters


% Compute the torque applied to the pendulum by gravity at each timestep
% Make sure to exclude the last theta value to match the length of omega_trimmed
tau = mass * g * length_to_COM * sin(theta_data(1:end-1));



disp("2.1.2");
% Define the fixed time-step (Delta t)
delta_t = 0.001; % Set the time-step to the value you have been given
% Use the diff function to compute the difference between each adjacent entry
omega_diff = diff(omega_data);
% Now compute y, which is the vector of approximate angular accelerations
y = omega_diff / delta_t;



disp("2.2.1");

% Time-step information
delta_t = 0.001; % seconds
t_final = 10.0; % seconds
t_grid = 0:delta_t:t_final; % time grid
N = numel(t_grid); % number of time steps



% Ensure dimensions match by trimming the last entry of omega_data
omega_trimmed = omega_data(1:end-1);

% Create the matrix A for the linear system, ensuring the dimensions match
A = [tau, omega_trimmed];

% Solve the linear system A * [c; d] = y
coefficients = A \ yk;

% Extract the estimated c and d from the solution
c = coefficients(1);
d = coefficients(2);

% Use c to calculate I (moment of inertia)
I = 1 / c;


disp("2.2.2");
% Define the fixed time-step (Delta t)
% Ensure theta_data and omega_data are column vectors
theta_data = theta_data(:);
omega_data = -omega_data(:);

% Compute the torque applied to the pendulum by gravity at each timestep
tau = mass * g * length_to_COM * sin(theta_data(1:end-1));

% Compute yk (angular acceleration), which has one fewer entry than omega_data
yk = diff(omega_data) / delta_t;

% Ensure dimensions match by trimming the last entry of omega_data to match yk
omega_trimmed = -omega_data(1:end-1);

% Create the matrix A for the linear system, ensuring the dimensions match
A = [tau, omega_trimmed];

% Solve the linear system A * [c; d] = yk
coefficients = A \ yk;

% Extract the estimated c and d from the solution
c_star = coefficients(1);
d_star = coefficients(2);
disp("C_star:")
disp(c_star);
disp("D_star:")
disp(d_star);


disp("2.2.3");
% Use c_star to calculate I (moment of inertia)
I_star = 1 / c_star;

% Output the results
fprintf('Estimated moment of inertia I: %f kg*m^2\n', I_star);
fprintf('Estimated damping coefficient d: %f N*m*s/rad\n', d_star);



disp("2.3");
c = 1 / I_star; % Constant based on known I


% Simulate or measure theta and omega data for the pendulum at the new location
% Here, we would use pendulum_far_from_home.simulateAndMeasure or similar method
% For example:
% [theta_data, omega_data] = pendulum_far_from_home.simulateAndMeasure(t_grid, initial_conditions);

% Calculate the differences between successive angular velocities
y = diff(omega_data) / delta_t;

tau_g = mass * length_to_COM * sin(theta_data(1:end-1)); % Note: Ensure theta_data is from the new simulation
A = [tau_g, -omega_data(1:end-1)];

% Solve the least squares problem
coefficients = A \ y;

% Extract g and d
g_star = coefficients(1) / (I_star); % Adjust for the known I_star
d_star = coefficients(2);

% Display estimated g and d
fprintf('Estimated g: %f m/s^2\n', g_star);
fprintf('Estimated d: %f N*m*s/rad\n', d_star);
