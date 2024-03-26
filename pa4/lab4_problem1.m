disp("1.2.3");
% Define the data points
x = [0; 1; 2; 3; 4];
y = [3.5; 4; 6; 5; 6.5];

% Create the plot
figure(1); clf;
hold on;
box on;
plot(x, y, '*');

% Define the slope (m_sol) and y-intercept (b_sol) computed previously
b_sol = 3.6;  % replace with the exact value you calculated
m_sol = 0.7;  % replace with the exact value you calculated

% Calculate the model's predictions
y_pred = m_sol * x + b_sol;

% Plot the model's predictions
plot(x, y_pred);

% Beautify the plot
title('Linear regression');
xlabel('x');
ylabel('y');
% Constrain the plot display to y > 0
ylim([0, inf]);



disp("1.3.1");
% Define the data points
x = [0; 1; 2; 3; 4];
y = [3.5; 4; 6; 5; 6.5];

% Create matrix A with a column of ones and a column with x values
A = [ones(size(x)), x];

% Define the least squares solution vector c* you computed previously
b_star = 3.6;  % replace with the exact value you calculated
m_star = 0.7;  % replace with the exact value you calculated
c_star = [b_star; m_star];

% Compute the predicted values A*c*
y_pred = A * c_star;

% Compute the residuals r = y - A*c*
r = y - y_pred;

% Verify that A*c* is orthogonal to r by checking if their dot product is zero
orthogonal = dot(A * c_star, r);

% Display the result
disp(['Dot product of A*c* and r: ', num2str(orthogonal)]);
if abs(orthogonal) < 1e-10  % Adjust the tolerance as necessary
    disp('A*c* is orthogonal to r.');
else
    disp('A*c* is not orthogonal to r.');
end




disp("1.3.2");
% Create matrix A with a column of ones and a column with x values
A = [ones(size(x)), x];


% Compute the residual vector r* = y - A*c*
r_star = y - A * c_star;

% Compute the norm of the residual vector r*
norm_r_star = norm(r_star);

% Display the norm of r*
disp(['Norm of r*: ', num2str(norm_r_star)]);

% Define visually estimated values for m and b
b_visual = 3.5; % replace with the value you estimated visually
m_visual = 0.75; % replace with the value you estimated visually
c_visual = [b_visual; m_visual];

% Compute the residual vector for visually estimated values of c
r_visual = y - A * c_visual;

% Compute the norm of the residual vector for visually estimated values of c
norm_r_visual = norm(r_visual);

% Display the norm of r for visually estimated c
disp(['Norm of r for visually estimated c: ', num2str(norm_r_visual)]);

% Compare the two norms
disp(['Difference in norms: ', num2str(norm_r_visual - norm_r_star)]);



disp("1.4.1");
% Define the validation data points
x_val = [0; 2; 3; 4; 5];
y_val = [3.1; 4.2; 4.9; 6.2; 6.9];


% Plot the original data and the validation data in a new figure
figure(2); clf;
hold on; box on;
plot(x, y, '*', 'DisplayName', 'Original Data');
plot(x_val, y_val, '*', 'DisplayName', 'Validation Data');

% Add labels and title
xlabel('x');
ylabel('y');
title('Model Testing with Validation Data');

% Add a legend to distinguish data sets
legend show;

% Constrain the plot display to y > 0
ylim([0, inf]);

% Optionally, plot the regression line using the model
x_combined = [x; x_val];
x_range = [min(x_combined), max(x_combined)];
y_model = m_star * x_range + b_star;
plot(x_range, y_model, 'DisplayName', 'Regression Line');

% Hold off the plot
hold off;




disp("1.4.2");
% Define the slope (m_star) and y-intercept (b_star) from your model
b_star = 3.6; % replace with the exact value you calculated
m_star = 0.7; % replace with the exact value you calculated

% Use the model to predict the y values for the validation data
y_val_pred = m_star * x_val + b_star;

% Plot the original data, validation data, and the linear regression predictions
figure(2); clf;
hold on; box on;
plot(x, y, '*', 'DisplayName', 'Training Data');
plot(x_val, y_val, '*', 'DisplayName', 'Validation Data');
plot(x_val, y_val_pred, '-', 'DisplayName', 'Linear Model (Predictions)');

% Beautify the plot
title('Linear regression (validation)');
xlabel('x');
ylabel('y');

% Set the legend for the three plots we added
lgd = legend('Training Data', 'Validation Data', 'Linear Model (Predictions)');
% Move the legend location
lgd.Location = 'northwest';

% Constrain the plot display to y > 0
ylim([0, inf]);

% Hold off the plot
hold off;




disp("1.4.3");
% Define the validation data points
x_val = [0; 2; 3; 4; 5];
y_val = [3.1; 4.2; 4.9; 6.2; 6.9];

% Create the validation matrix A_val with a column of ones and a column with x_val values
A_val = [ones(size(x_val)), x_val];

% Use the least squares solution vector c* you computed previously
b_star = 3.6; % replace with the exact value you calculated
m_star = 0.7; % replace with the exact value you calculated
c_star = [b_star; m_star];

% Compute the residual vector r_val = y_val - A_val*c*
r_val = y_val - A_val * c_star;

% Compute the norm of the residual vector r_val
norm_r_val = norm(r_val);

% Display the norm of r_val
disp(['Norm of r_val: ', num2str(norm_r_val)]);

% You should also have the norm of r_star from the training data computed previously
% norm_r_star = ...

% Display the comparison
% disp(['Difference in norms: ', num2str(norm_r_val - norm_r_star)]);

% Depending on your results, provide commentary on model generalization
if norm_r_val > norm_r_star
    disp('The validation error is larger than the training error. The model may not generalize as well as desired.');
elseif norm_r_val < norm_r_star
    disp('The validation error is smaller than the training error. The model generalizes well.');
else
    disp('The validation error is equal to the training error. The model generalizes consistently.');
end



disp("1.5.1");
% Define the training data points
x = [0; 1; 2; 3; 4];
y = [3.5; 4; 6; 5; 6.5];

% Create the Vandermonde matrix X for the polynomial fitting
X = fliplr(vander(x));

% Solve for the polynomial coefficients a
a = X \ y; % This uses MATLAB's backslash operator to solve the linear system

% Display the polynomial coefficients
disp('The polynomial coefficients a are:');
disp(a);

% Create a function handle for the resulting polynomial interpolant p(x)
p = @(x) a(1) + a(2)*x + a(3)*x.^2 + a(4)*x.^3 + a(5)*x.^4;

% Display the polynomial interpolant
disp('The resulting polynomial interpolant p(x) is:');
disp(['p(x) = ', num2str(a(1)), ' + ', num2str(a(2)), 'x + ', ...
      num2str(a(3)), 'x^2 + ', num2str(a(4)), 'x^3 + ', num2str(a(5)), 'x^4']);



disp("1.5.2");

% Create a grid of x coordinates to use for plotting the polynomial
x_grid = linspace(0, 5);

% Evaluate p(x) at all the points in x_grid
p_of_x_grid = a(1) + a(2) * x_grid + a(3) * x_grid.^2 + a(4) * x_grid.^3 + a(5) * x_grid.^4;

% Define the validation data points
x_val = [0; 2; 3; 4; 5];
y_val = [3.1; 4.2; 4.9; 6.2; 6.9];

% Define the slope (m_star) and y-intercept (b_star) from your linear model
% Note: Replace these with the values you found previously
m_star = 0.7;
b_star = 3.6;

% Predictions using the linear model
y_val_pred = m_star * x_val + b_star;

% Create a new figure and plot the training data, validation data, and the predictions
figure(3); clf;
hold on;
box on;

% Plot given training data
plot(x, y, '*', 'DisplayName', 'Training Data');

% Plot validation data
plot(x_val, y_val, '*', 'DisplayName', 'Validation Data');

% Plot our least-squares predictions
plot(x_val, y_val_pred, 'DisplayName', 'Linear Model (Predictions)');

% Plot the polynomial interpolation
plot(x_grid, p_of_x_grid, 'DisplayName', 'Polynomial Model (Predictions)');

% Beautify the plot
xlabel('x');
ylabel('y');
title('Polynomial Interpolation vs. Linear Least Squares');

% Add a legend to label each plot
lgd = legend('Training Data', 'Validation Data', 'Linear Model (Predictions)', 'Polynomial Model (Predictions)');
lgd.Location = 'northwest';

% Adjust the limits
ylim([0, 10]);

% Display the plot
hold off;

