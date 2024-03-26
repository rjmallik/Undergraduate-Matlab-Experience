disp("Problem 3.1");
%Part a all on lab report
disp("Part B:");
%Part b
% Given values
R1 = 100; % Ohms
R2 = 150; % Ohms
R3 = 200; % Ohms
R4 = 50;  % Ohms
V1 = 1.5; % Volts
V2 = 1.5; % Volts

% Define the coefficient matrix A and the vector b
% The matrix A should be based on the chosen direction of the loop currents and the given equations.
A = [1, 0, 0, 1, 0, 0;     % i1 = i2 + i4 (from equation 3)
     -1, 0, 0, 1, 0, 1;    % i4 + i6 = i1 (from equation 4)
     R1, 0, 0, R4, 0, 0;   % V1 - R1*i1 - R4*i4 = 0 (from equation 5)
     0, R2, 0, -R4, 0, 0;  % V2 + R4*i4 - R2*i2 = 0 (from equation 6, note the sign change for i4)
     0, 1, 0, 0, -1, 0;    % i2 = i3 + i5 (from equation 8)
     0, 0, 1, 0, 1, -1];   % i5 = i6 + i3 (from equation 9)

b = [0; 0; V1; V2; 0; 0]; % Corresponding values on the right side of the equations

% Solve the system using LU decomposition
[L, U, P] = lu(A);  % LU decomposition with permutation
y = L \ (P * b);    % Forward substitution
i_vector = U \ y;    % Backward substitution

% Display the currents
disp('Currents are:');
disp(['i1 = ' num2str(i_vector(1)) ' A']);
disp(['i2 = ' num2str(i_vector(2)) ' A']);
disp(['i3 = ' num2str(i_vector(3)) ' A']);
disp(['i4 = ' num2str(i_vector(4)) ' A']);
disp(['i5 = ' num2str(i_vector(5)) ' A']);
disp(['i6 = ' num2str(i_vector(6)) ' A']);

% Verify the solution
residual = A * i_vector - b;
disp('Residuals should be close to zero to confirm the solution is accurate:');
disp(residual);

disp("Part C:")
% Part c(checked to make sure the values match up with my equations)
% Update the voltage source value for V2
% Define the given values
R1 = 100; % Ohms
R2 = 150; % Ohms
R3 = 200; % Ohms
R4 = 50;  % Ohms
V1 = 1.5; % Volts for part (b)
V2 = 9;   % Volts for part (c)

% Define the coefficient matrix A and vector b for part (c)
A = [1, -1, 0, -1, 0, 0;   % i1 = i2 + i4
     0, 1, -1, 0, -1, 0;   % i2 = i3 + i5
     0, 0, 1, 0, 1, -1;    % i5 = i6 + i3
     R1, 0, 0, R4, 0, 0;   % V1 - R1*i1 - R4*i4 = 0
     0, R2, 0, -R4, 0, 0;  % V2 - R2*i2 + R4*i4 = 0
     0, R2, R3, 0, 0, 0];  % -V2 + R2*i2 + R3*i3 = 0

b = [0; 0; 0; V1; V2; -V2]; % Update the vector b for the new voltage V2 for part (c)

% Solve the system using LU decomposition
[L, U, P] = lu(A); % LU decomposition with permutation
y = L \ (P * b);   % Forward substitution using the new b
i_vector_c = U \ y; % Backward substitution to solve for the currents

% Display the new currents for part (c)
disp('New currents with V2 = 9V are:');
disp(i_vector_c);


disp("Problem 3.2");
disp("Part A:")
% Part a
N = 200; % Number of points
t = linspace(0, 20, N); % Create a row vector of time points from 0 to 20
disp(t);
% Initialize voltage arrays
V1 = sin(t); % V1(t) as the sine of t
V2 = cos(t); % V2(t) as the cosine of t
% Parb b
V = zeros(N, 2);
i = zeros(N, 6);
% Part c & d
% Given resistor values from the previous problem
R1 = 100; % Ohms
R2 = 150; % Ohms
R3 = 200; % Ohms
R4 = 50;  % Ohms

% Define the matrix A using the resistor values
A = [1, -1, 0, -1, 0, 0;   % i1 = i2 + i4
     0, 1, -1, 0, -1, 0;   % i2 = i3 + i5
     0, 0, 1, 0, 1, -1;    % i5 = i6 + i3
     R1, 0, 0, R4, 0, 0;   % V1 - R1*i1 - R4*i4 = 0
     0, R2, 0, -R4, 0, 0;  % V2 - R2*i2 + R4*i4 = 0
     0, R2, R3, 0, 0, 0];  % -V2 + R2*i2 + R3*i3 = 0

% Precompute the LU factorization of A for efficiency
[L, U, P] = lu(A);

% Initialize N, t, V, and i arrays
N = 200;
t = linspace(0, 20, N);
V = zeros(N, 2); % To store V1(t) and V2(t) at each time point
i = zeros(N, 6); % To store i1(t) to i6(t) at each time point

% Loop through each time value
for k = 1:N
    % Compute V1(t) and V2(t)
    V(k, 1) = sin(t(k));
    V(k, 2) = cos(t(k));
    
    % Redefine b for the current time step
    b = [0; 0; 0; V(k, 1); V(k, 2); -V(k, 2)];
    
    % Solve the system for i(t) using the precomputed LU factors
    y = L \ (P * b);   % Forward substitution
    i(k, :) = U \ y;   % Backward substitution to find i at time t(k)
end

% At this point, V contains all voltage values and i contains all current values for each time t

% Part e
% After computing V and i arrays
figure (1) ; clf
% Plot the voltages vs. time
subplot(2, 1, 1); % This creates a subplot in the first row
plot(t, V); % Plotting both V1(t) and V2(t) vs. time
title('Voltages vs. Time');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V1(t)', 'V2(t)'); % Add a legend to distinguish between V1 and V2

% Plot the currents vs. time
subplot(2, 1, 2); % This creates a subplot in the second row
plot(t, i); % Plotting all currents i1(t) to i6(t) vs. time
title('Currents vs. Time');
xlabel('Time (s)');
ylabel('Current (A)');
legend('i1(t)', 'i2(t)', 'i3(t)', 'i4(t)', 'i5(t)', 'i6(t)'); % Add a legend to distinguish between currents

% Adjust layout to make sure everything fits nicely
sgtitle('Circuit Analysis Over Time'); % Super title for the whole figure

